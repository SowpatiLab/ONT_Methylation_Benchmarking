rule p5_to_b5:
    input: "datasets/pod5/{experiment}_{sr}kHz"
    output: directory("datasets/blow5/{experiment}_{sr}kHz")
    threads: 20
    log: "log/pod5_to_blow5/{experiment}_{sr}kHz.log"
    shell: "mkdir -p {output}; {TIME} blue-crab p2s {input} -t {threads} -o {output}/{wildcards.experiment}_{wildcards.sr}kHz.blow5 2> {log}"

rule f5c_idx_fastq:
    input:  
        blow5="datasets/blow5/{experiment}_{sr}kHz",
        fastq="fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq"
    output: 
        fastq=multiext("fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq", '.index', '.index.fai', '.index.gzi')
    threads: 20
    log: "log/f5c_idx/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell:  "{TIME} /data1/tools/f5c-v1.5/f5c_x86_64_linux_cuda index --slow5 {input.blow5}/{wildcards.experiment}_{wildcards.sr}kHz.blow5 {input.fastq}; 2> {log}"

rule f5c_call:
    input: 
        blow5="datasets/blow5/{experiment}_{sr}kHz",
        bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
        fastq="fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq",
        fastq_idx="fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq.index"
    output: "f5c/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    threads: 20
    priority: 1
    resources: gpu=1
    params: 
        ref=getRef
    log: "log/f5c_call/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: "{log} /data1/tools/f5c-v1.5/f5c_x86_64_linux_cuda call-methylation \
                --cuda-dev-id 1 \
                -x hpc-high \
                -K 2048 \
                -B 40.0M \
                -t {threads} \
                --slow5 {input.blow5}/{wildcards.experiment}_{wildcards.sr}kHz.blow5 \
                -b {input.bam} \
                -g {params.ref}.fa \
                -r {input.fastq} \
                --pore r10 > {output}"

rule f5c_restrand:
    input:  'f5c/readwise/{experiment}_5kHz_{acc}_v{ver}.tsv'
    output: 'f5c/readwise_stranded/{experiment}_5kHz_{acc}_v{ver}.tsv'
    log:    'log/f5c_restrand/{experiment}_5kHz_{acc}_v{ver}.tsv'
    threads: 20
    params: ref=getRef
    shell: "{TIME} python workflow/scripts/f5c_restrand.py -i {input} -r {params.ref}.fa -o {output} 2> {log}"

rule f5c_aggregate:
    input:  "f5c/readwise{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    output: "f5c/aggregated{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    log: "log/f5c_aggregate{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    threads: 20
    shell: "{TIME} /data1/tools/f5c-v1.5/f5c_x86_64_linux_cuda  meth-freq -i {input} -s > {output} 2> {log}"

rule f5c_fasta:
    input:  "f5c/aggregated{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    output: "f5c/aggregated_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    log:    "log/f5c_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        ref=getRef
    threads: 20
    run: 
        if wildcards.strandstate=="_stranded":
            shell('''
                tmp=$(mktemp /tmp/f5c_bed_fasta.XXXX);
                
                bedtools getfasta -fi {params.ref}.fa -bed <(tail -n +2 {input}) -bedOut | \
                    awk 'BEGIN{{OFS="\t"}} {{ if(toupper($9)=="C") $9="+"; else if(toupper($9)=="G") $9="-"; print $1,$2,$3,$5,".",$9,$6,$5-$6,$7 }}' > $tmp;

                {TIME} bedtools slop -l 5 -r 11 -s -g {params.ref}.genome -i $tmp \
                    | bedtools getfasta -fi {params.ref}.fa -bed - -tab -s | cut -f 2 | paste $tmp - > {output} 2> {log}
            ''')
        else:
            shell('''
                tmp=$(mktemp /tmp/f5c_bed_fasta.XXXX);

                sed '1d' {input} | awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$2+1,$5,".","+",$6,$5-$6,$7}}' > $tmp;
                {TIME} bedtools slop -l 5 -r 11 -s -g {params.ref}.genome -i $tmp \
                    | bedtools getfasta -fi {params.ref}.fa -bed - -tab -s | cut -f 2 | paste $tmp - > {output} 2> {log}
            ''')
            
rule consolidate_f5c:
    input:  "f5c/aggregated_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    output: "meta/f5c{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    threads: 20
    run:
        import polars as pl
        t = pl.read_csv(input[0], separator='\t', has_header=False, new_columns=['chrom', 'p1', 'p2', 'coverage', 's', 'strand', 'M', 'UM', 'per', 'full_context'])
        select  = ['chrom', 'p1', 'p2', 'mod', 'coverage', 'strand', 'M', 'UM', 'per', 'flowcell', 'tool', 'model', 'sample', 'acc', 'sample_rate', 'species', 'full_context']
        species = wildcards.experiment.split('_')[0]
        
        t = (
            t.with_columns([
                pl.lit('5mC').alias('mod'),
                pl.lit('r10.4.1').alias('flowcell'),
                (pl.col('M')+pl.col('UM')).alias('coverage'),
                pl.lit(wildcards.experiment).alias('sample'),
                pl.lit(f'f5c{wildcards.strandstate}').alias('tool'),
                pl.lit(wildcards.acc).alias('acc'),
                pl.lit('.').alias('model'),
                pl.lit(f"{wildcards.sr}kHz").alias('sample_rate'),
                pl.lit(species).alias('species'),
                (pl.col('per')*100).alias('per')
            ])
        )
        t.select(select).write_csv(output[0], separator='\t', include_header=True)