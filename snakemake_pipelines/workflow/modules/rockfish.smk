rule rockfish:
    input: 
        pod5="datasets/pod5/{experiment}_{sr}kHz",
        bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: 
        main=directory("rockfish/{experiment}_{sr}kHz_{acc}_v{ver}"),
        pred="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}/predictions.tsv",
        # ref_map="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv"
    threads: 40
    priority: 1
    resources: gpu=1
    log: "log/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    conda: 'rockfish_flash'
    params: 
        model="/data1/tools/rockfishmodels/rf_5kHz.ckpt",
        ref=getRef
    shell: '''
        {log} mkdir -p {output.main}; 
        cd {output.main}; 
        /usr/bin/time -v -o ../../{log} rockfish inference \
            -i ../../{input.pod5} \
            --bam_path ../../{input.bam} \
            --model_path {params.model} \
            -r -t {threads} \
            -b 16384 -d 1;
        cd ../../;
        '''

rule rockfish_map:
    input:  "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv"
    threads: 40
    priority: 1
    resources: gpu=1
    log: "log/rockfish_map/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    conda: 'rockfish_flash'
    params:
        ref=getRef
    shell: "{TIME} python /data1/tools/rockfish/scripts/extract_ref_pos.py --workers {threads} {input} {params.ref}.fa > {output} 2> {log}"

rule rockfish_intersect:
    input: 
        ref="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv",
        pre="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}/predictions.tsv"
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_remapped.tsv"
    threads: 20
    log: "log/rockfish_intersect/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    priority: 1
    shell: "{TIME} python workflow/scripts/rockfish_intersect.py -i {input.pre} -r {input.ref} -o {output} 2> {log}"

rule rockfish_aggregate:
    input:  "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_remapped.tsv"
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated.tsv"
    threads: 20
    log: "log/rockfish_aggregate/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: "{TIME} python workflow/scripts/rockfish_aggregate.py -i {input} -o {output} 2> {log}"

rule rockfish_getfasta:
    input:  "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated.tsv"
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated_reference.tsv"
    threads: 20
    log: "log/rockfish_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: ref=getRef
    shell: '''
        tmp=$(mktemp /tmp/rockfish_bed_fasta.XXXX);
        bedtools getfasta -fi {params.ref}.fa -bed {input} -bedOut | \
            awk 'BEGIN{{OFS="\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > $tmp;
        {TIME} bedtools slop -g {params.ref}.genome -l 5 -r 11 -s -i $tmp | bedtools getfasta -s -fi {params.ref}.fa -bed - -tab | cut -f 2 | paste $tmp - > {output} 2> {log}
    '''

rule consolidate_rockfish:
    input:  "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated_reference.tsv"
    output: "meta/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.tsv"
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
                pl.lit(wildcards.experiment).alias('sample'),
                pl.lit('rockfish').alias('tool'),
                pl.lit(wildcards.acc).alias('acc'),
                pl.lit('.').alias('model'),
                pl.lit(f"{wildcards.sr}kHz").alias('sample_rate'),
                pl.lit(species).alias('species'),
                (pl.col('per')*100).alias('per')
            ])
        )
        t.select(select).write_csv(output[0], separator='\t', include_header=True)
