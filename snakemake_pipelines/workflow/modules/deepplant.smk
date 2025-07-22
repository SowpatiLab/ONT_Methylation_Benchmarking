if is_bacterial==False:
    rule deep_plant_split:
        input:  
            pod5="p5_splits/signal/{experiment}_{sr}kHz/{p5}",
            bam="bam/sorted_move_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}.bam",
        output: 
            main=directory("deepplant/readwise_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}"),
            cpg="deepplant/readwise_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}/cpg_result.txt",
            chg="deepplant/readwise_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}/chg_result.txt",
            chh="deepplant/readwise_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}/chh_result.txt"
        threads: 32
        log: "log/deeplant_split/{experiment}_{sr}kHz_{acc}_v{ver}_v{ver}_{p5}.log"
        conda: "deepplant"
        priority: 2
        params: 
            ref=getRef,
            model_dir="/data1/tools/DeepPlant/model/bilstm/"
        shell: '''
            mkdir {output} -p;
            /data1/tools/DeepPlant/build/DeepPlant extract_and_call_mods \
                    {input.pod5} {input.bam} {params.ref}.fa DNA {output.main} {params.model_dir} \
                    51 \
                    51 \
                    13 \
                    8 \
                    4 \
                    2048 
        '''

    rule deepplant_cat_readwise:
        input: 
            cpg=expand("deepplant/readwise_split/{{experiment}}_{{sr}}kHz_{{acc}}_v{{ver}}/{p5}/cpg_result.txt", p5=P5S),
            chg=expand("deepplant/readwise_split/{{experiment}}_{{sr}}kHz_{{acc}}_v{{ver}}/{p5}/chg_result.txt", p5=P5S),
            chh=expand("deepplant/readwise_split/{{experiment}}_{{sr}}kHz_{{acc}}_v{{ver}}/{p5}/chh_result.txt", p5=P5S)
        output: 
            main=directory("deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}"),
            cpg="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/cpg_result.txt",
            chg="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/chg_result.txt",
            chh="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/chh_result.txt",
        log: "log/deepplant_concat/{experiment}_{sr}kHz_{acc}_v{ver}"
        threads: 20
        shell:'''
            [ ! -d {output.main} ] && mkdir -p {output.main}
            ( head -n1 -q {input.cpg} | head -n 1; tail -n +2 -q {input.cpg} ) > {output.cpg};
            ( head -n1 -q {input.chg} | head -n 1; tail -n +2 -q {input.chg} ) > {output.chg};
            ( head -n1 -q {input.chh} | head -n 1; tail -n +2 -q {input.chh} ) > {output.chh};
        '''
else:
    rule deep_plant:
        input:  
            pod5="datasets/pod5/{experiment}_{sr}kHz",
            bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
        output: 
            main=directory("deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}"),
            cpg="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/cpg_result.txt",
            chg="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/chg_result.txt",
            chh="deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/chh_result.txt"
        threads: 32
        log: "log/deeplant/{experiment}_{sr}kHz_{acc}_v{ver}_v{ver}.log"
        conda: "deepplant"
        priority: 3
        params: 
            ref=getRef,
            model_dir="/data1/tools/DeepPlant/model/bilstm/"
        shell: '''
            mkdir {output} -p;
            {TIME} /usr/bin/time -v -o {log} /data1/tools/DeepPlant/build/DeepPlant extract_and_call_mods \
                    {input.pod5} {input.bam} {params.ref}.fa DNA {output.main} {params.model_dir} \
                    51 \
                    51 \
                    13 \
                    8 \
                    4 \
                    2048 2> {log}
        '''

rule deepplant_aggregate:
    input:  "deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}/{context}_result.txt"
    output: "deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    threads: 20
    log:    "log/deepplant_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell:  "{TIME}   python workflow/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} \
                --threshold 0.5 2> {log}"

rule deepplant_rebed:
    input:  "deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    output: "deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    threads: 20
    log:    "log/deepplant_rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: ref=getRef
    shell: '''
        {TIME} bedtools getfasta -fi {params.ref}.fa -bed {input} -bedOut | \
            awk 'BEGIN{{OFS="\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output} 2> {log};
    '''

rule deepplant_addfasta:
    input:  "deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    output: "deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.out"
    threads: 20
    params: ref=getRef
    log:    "log/deepplant_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell: "{TIME} bedtools slop -g {params.ref}.genome -l 5 -r 11 -i {input} -s | bedtools getfasta -fi {params.ref}.fa -bed - -tab -s | cut -f 2 | paste {input} - > {output} 2> {log}"

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

