
rule setup_f5c:
    output: "{tooling}/setup_status/f5c_{version}.status"
    params: version="v1.5"
    shell: """
        [ ! -d {wildcards.tooling} ] && mkdir {wildcards.tooling} -p
        cd {wildcards.tooling}
        
        wget "https://github.com/hasindu2008/f5c/releases/download/{params.version}/f5c-{params.version}-binaries.tar.gz" && \
            tar xvf f5c-{params.version}-binaries.tar.gz
        
        cd ../
        echo "f5c dir: $(realpath {wildcards.tooling}/f5c-{params.version})" > {output}
    """

rule p5_to_b5:
    input: config['pod5_dir'] + "/{experiment}_{sr}kHz"
    output: directory("datasets/blow5/{experiment}_{sr}kHz")
    threads: 20
    log: "log/pod5_to_blow5/{experiment}_{sr}kHz.log"
    shell: ntsh("mkdir -p {output}; {TIME} blue-crab p2s {input} -t {threads} -o {output}/{wildcards.experiment}_{wildcards.sr}kHz.blow5")

rule f5c_idx_fastq:
    input:  
        # setup=expand('{tooling}/setup_status/f5c_{version}.status', tooling=config['tooling_dir'], version=config['tooling']['f5c']['version']),
        blow5="datasets/blow5/{experiment}_{sr}kHz",
        fastq="fastq_nanoq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq"
    output: 
        fastq=multiext("fastq_nanoq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq", '.index', '.index.fai', '.index.gzi')
    threads: 20
    log: "log/f5c_idx/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params:
        execut=lambda x: f"{config['tooling_dir']}/{config['tooling']['f5c']['install']}/{config['tooling']['f5c']['execut']}"
    shell:  sh("{params.execut} index --slow5 {input.blow5}/{wildcards.experiment}_{wildcards.sr}kHz.blow5 {input.fastq};")

rule f5c_call:
    input: 
        # setup=expand('{tooling}/setup_status/f5c_{version}.status', tooling=config['tooling_dir'], version=config['tooling']['f5c']['version']),
        blow5="datasets/blow5/{experiment}_{sr}kHz",
        bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
        fastq="fastq_nanoq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq",
        fastq_idx="fastq_nanoq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq.index"
    output: "f5c/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    threads: 20
    priority: 1
    resources: gpu=1
    params: 
        ref=getRef,
        execut=lambda x: f"{config['tooling_dir']}/{config['tooling']['f5c']['install']}/{config['tooling']['f5c']['execut']}"
    log: "log/f5c_call/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("{params.execut} call-methylation \
                --cuda-dev-id 1 \
                -x hpc-high \
                -K 2048 \
                -B 40.0M \
                -t {threads} \
                --slow5 {input.blow5}/{wildcards.experiment}_{wildcards.sr}kHz.blow5 \
                -b {input.bam} \
                -g {params.ref} \
                -r {input.fastq} \
                --pore r10 > {output}")

rule f5c_restrand:
    input:  'f5c/readwise/{experiment}_5kHz_{acc}_v{ver}.tsv'
    output: 'f5c/readwise_stranded/{experiment}_5kHz_{acc}_v{ver}.tsv'
    log:    'log/f5c_restrand/{experiment}_5kHz_{acc}_v{ver}.tsv'
    threads: 20
    params: ref=getRef
    shell: sh("python workflow/scripts/f5c_restrand.py -i {input} -r {params.ref} -o {output}")

rule f5c_aggregate:
    input:  
        # setup=expand('{tooling}/setup_status/f5c_{version}.status', tooling=config['tooling_dir'], version=config['tooling']['f5c']['version']),
        out="f5c/readwise{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    output: "f5c/aggregated{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    log: "log/f5c_aggregate{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    threads: 20
    params:
        execut=lambda x: f"{config['tooling_dir']}/{config['tooling']['f5c']['install']}/{config['tooling']['f5c']['execut']}"
    shell: sh("{params.execut} meth-freq -i {input.out} -s > {output}")

rule f5c_fasta:
    input:  
        bed="f5c/aggregated{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: "f5c/aggregated_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    log:    "log/f5c_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    threads: 20
    run: 
        if wildcards.strandstate=="_stranded":
            shell(sh('''
                tmp=$(mktemp /tmp/f5c_bed_fasta.XXXX);
                
                {BEDTOOLS} getfasta -fi {input.fasta} -bed <(tail -n +2 {input.bed}) -bedOut | \
                    awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($9)=="C") $9="+"; else if(toupper($9)=="G") $9="-"; print $1,$2,$3,$5,".",$9,$6,$5-$6,$7 }}' > $tmp;

                {TIME} {BEDTOOLS} slop -l 5 -r 11 -s -g {input.genome} -i $tmp \
                    | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $tmp - > {output};
                rm $tmp
            '''))
        else:
            shell(ntsh('''
                tmp=$(mktemp /tmp/f5c_bed_fasta.XXXX);

                sed '1d' {input.bed} | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$2+1,$5,".","+",$6,$5-$6,$7}}' > $tmp;
                {TIME} {BEDTOOLS} slop -l 5 -r 11 -s -g {input.genome} -i $tmp \
                    | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $tmp - > {output};
                rm $tmp
            '''))
            
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