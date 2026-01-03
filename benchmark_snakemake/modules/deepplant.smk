rule bam_fn_reorder:
    input:  "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    output: "bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    log:    "log/sorted_move_by_fntag/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    conda: "snakemake"
    threads: 20
    resources:
        queue='hm-q',
        mem='250gb'
    shell:  sh("{SAMTOOLS} sort -t fn {input} -@ {threads} -O BAM -o {output}")

rule deepplant_call:
    input:  
        pod5=config['pod5_dir'] + "/{experiment}_{sr}kHz",
        bam="bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: 
        main=directory("deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}"),
        cpg="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/cpg_result.txt",
        chg="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/chg_result.txt",
        chh="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/chh_result.txt"
    threads: 32
    log: "log/deeplant/{experiment}_{sr}kHz_{acc}_v{ver}_v{ver}.log"
    conda: config['tooling']['deepplant']['conda']
    priority: 3
    params: 
        ref=getRef,
        model_dir=config['tooling']['deepplant']['model'],
        execut=config['tooling']['deepplant']['execut']
    shell: ntsh('''
        mkdir {output.main} -p;
        {TIME} {params.execut} extract_and_call_mods \
            {input.pod5} {input.bam} {params.ref} DNA {output.main} {params.model_dir} \
            51 51 13 \
            8 4 2048
    ''')
    # shell: '''
    #     mkdir {output} -p;
    #     {TIME} singularity run --bind $(realpath ./):/mnt/work/ --nv /data1/ccmb/reuben/benchmarking/deeptools.sif deepplant \
    #             {input.pod5} {input.bam} {params.ref} DNA /mnt/work/{output.main} /deepplant/model/bilstm/ \
    #             51 \
    #             51 \
    #             13 \
    #             8 \
    #             4 \
    #             2048 2> >(tee -a {log} >&2
    # '''

rule deepplant_aggregate:
    input:  "deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/{context}_result.txt"
    output: "deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    threads: 20
    # conda:  "snakemake"
    log:    "log/deepplant_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell:  sh("python workflow/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} \
                --threshold 0.5")

rule deepplant_rebed:
    input:  "deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    output: "deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    threads: 20
    log:    "log/deepplant_rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: ref=getRef
    shell: sh('''
        {BEDTOOLS} getfasta -fi {params.ref} -bed {input} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output};
    ''')

rule deepplant_addfasta:
    input: 
        bed="deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: "deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    threads: 20
    log:    "log/deepplant_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell: sh("{BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i {input.bed} -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste {input.bed} - > {output}")

rule consolidate_deepplant:
    input:  "deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    output: "meta/deepplant/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.tsv"
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
                pl.lit('deepplant').alias('tool'),
                pl.lit(wildcards.acc).alias('acc'),
                pl.lit('2024').alias('model'),
                pl.lit(f"{wildcards.sr}kHz").alias('sample_rate'),
                pl.lit(species).alias('species')
            ])
        )
        t.select(select).write_csv(output[0], separator='\t', include_header=True)

