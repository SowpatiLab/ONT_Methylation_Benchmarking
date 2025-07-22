
if is_bacterial==False:
    rule deep_bam_split:
        input:  
            pod5="p5_splits/signal/{experiment}_{sr}kHz/{p5}",
            bam="bam/sorted_move_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}.bam",
        output: "deepbam/readwise_split/{experiment}_{sr}kHz_{acc}_v{ver}/{p5}.out"
        threads: 32
        priority: 0
        resources: gpu=1
        conda: 'deepbam'
        log: "log/deepbam_split/{experiment}_{sr}kHz_{acc}_v{ver}_{p5}.log"
        params: 
            ref=getRef,
            model="/home/tej/zips/DeepBAM/traced_script_module/LSTM_20240524_newfeature_script_b9_s15_epoch25_accuracy0.9742.pt"
        shell: '''{TIME} DeepBAM extract_and_call_mods \
            {input.pod5} {input.bam} {params.ref}.fa \
            DNA {output} {params.model} \
            51 \
            8 \
            4 \
            1024 \
            CG \
            0 2> {log};
        '''

    rule deepbam_cat_readwise:
        input: expand("deepbam/readwise_split/{{experiment}}_{{sr}}kHz_{{acc}}_v{{ver}}/{p5}.out", p5=P5S),
        output: 
            main=directory("deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}"),
            out="deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.out",
        log: "log/deepbam_concat/{experiment}_{sr}kHz_{acc}_v{ver}"
        threads: 32
        shell:'''
            [ ! -d {output.main} ] && mkdir -p {output.main}
            {TIME}  ( head -n1 -q {input} | head -n 1; tail -n +2 -q {input} ) > {output.out} 2> {log};
        '''
else: 
    rule deep_bam:
            input:  
                pod5="datasets/pod5/{experiment}_{sr}kHz/",
                bam="bam/sorted_move/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
            output: "deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.out"
            threads: 32
            priority: 0
            resources: gpu=1
            conda: 'deepbam'
            log: "log/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.log"
            params: 
                ref=getRef,
                # model="/home/tej/zips/DeepBAM/traced_script_module/LSTM_20231111_b51_s15_epoch4_acc9883.pt"
                model="/home/tej/zips/DeepBAM/traced_script_module/LSTM_20240524_newfeature_script_b9_s15_epoch25_accuracy0.9742.pt"
            shell: '''{TIME}  DeepBAM extract_and_call_mods \
                {input.pod5} {input.bam} {params.ref}.fa \
                DNA {output} {params.model} \
                51 \
                8 \
                4 \
                1024 \
                CG \
                0 2> {log};
            '''

rule deepbam_aggregate:
    input:  "deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    output: "deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 20
    log:    "log/deepbam_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell:  "{TIME}   python workflow/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} \
                --threshold 0.5 2> {log}"

rule deepbam_rebed:
    input:  "deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    output: "deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 20
    log:    "log/deepbam_rebed/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: ref=getRef
    shell: '''
        {TIME}   bedtools getfasta -fi {params.ref}.fa -bed {input} -bedOut | \
            awk 'BEGIN{{OFS="\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output} 2> {log};
    '''

rule deepbam_addfasta:
    input:  "deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    output: "deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.out"
    threads: 20
    params: ref=getRef
    log:    "log/deepbam_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: "{TIME} bedtools slop -g {params.ref}.genome -l 5 -r 11 -i {input} -s | bedtools getfasta -fi {params.ref}.fa -bed - -tab -s | cut -f 2 | paste {input} - > {output} 2> {log}"


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
