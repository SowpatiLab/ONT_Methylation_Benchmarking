rule deepbam_call:
    input:  
        pod5=config['pod5_dir'] + "/{experiment}_{sr}kHz/",
        bam="bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: "deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 32
    priority: 0
    resources: gpu=1
    conda: config['tooling']['deepbam']['conda']
    log: "log/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        ref=getRef,
        model=config['tooling']['deepbam']['model'],
        execut=config['tooling']['deepbam']['execut']
    shell: ntsh('''
        echo $CONDA_PREFIX;
        {TIME} {params.execut} extract_and_call_mods \
        {input.pod5} {input.bam} {params.ref} \
        DNA {output} {params.model} \
        51 \
        8 \
        4 \
        1024 \
        CG \
        0
    ''')

rule deepbam_aggregate:
    input:  "deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    output: "deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 20
    log:    "log/deepbam_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell:  sh("python workflow/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} \
                --threshold 0.5")

rule deepbam_rebed:
    input:  
        bed="deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.out",
        fasta=lambda wildcards: getRef(wildcards)
    output: "deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 20
    log:    "log/deepbam_rebed/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: ref=getRef
    shell: sh('''
        {BEDTOOLS} getfasta -fi {input.fasta} -bed {input.bed} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output};
    ''')

rule deepbam_addfasta:
    input:  
        bed="deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.out",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: "deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 20
    log:    "log/deepbam_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("{BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i {input.bed} -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste {input.bed} - > {output}")

rule consolidate_deepbam:
    input:  "deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    output: "meta/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
    threads: 20
    run:
        import polars as pl
        t = pl.read_csv(input[0], separator='\t', has_header=False, new_columns=['chrom', 'p1', 'p2', 'coverage', 's', 'strand', 'M', 'UM', 'per','full_context'])
        select  = ['chrom', 'p1', 'p2', 'mod', 'coverage', 'strand', 'M', 'UM', 'per', 'flowcell', 'tool', 'model', 'sample', 'acc', 'sample_rate', 'species', 'full_context']
        species = wildcards.experiment.split('_')[0]
        
        t = (
            t.with_columns([
                pl.lit('5mC').alias('mod'),
                pl.lit('r10.4.1').alias('flowcell'),
                (pl.col('M')+pl.col('UM')).alias('coverage'),
                pl.lit(wildcards.experiment).alias('sample'),
                pl.lit('deepbam').alias('tool'),
                pl.lit(wildcards.acc).alias('acc'),
                pl.lit('2024').alias('model'),
                pl.lit(f"{wildcards.sr}kHz").alias('sample_rate'),
                pl.lit(species).alias('species')
            ])
        )
        t.select(select).write_csv(output[0], separator='\t', include_header=True)
