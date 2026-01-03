rule setup_rockfish:
    output: "{tooling}/setup_status/rockfish.status"
    params: 
        conda_env=config['tooling']['rockfish']['conda'],
        model=subpath(config['tooling']['rockfish']['model'], parent=True)
    shell: """
        source "$(conda info --base)/etc/profile.d/conda.sh"
        [ ! -d {wildcards.tooling} ] && mkdir {wildcards.tooling} -p
        cd {wildcards.tooling}
        conda create --name {params.conda_env} python=3.9 -y
        conda activate {params.conda_env}
        git clone -b r10.4.1 https://github.com/lbcb-sci/rockfish.git --single-branch rockfish && cd rockfish
        pip install --extra-index-url https://download.pytorch.org/whl/cu118 .

        cd ../
        [ ! -d {params.model} ] && mkdir -p {params.model}
        rockfish download -m 5kHz -s {params.model}

        echo rockfish dir: $(realpath {wildcards.tooling}/rockfish) >  {output}
        echo rockfish model dir: {params.model} >> {output}
    """

rule rockfish_call:
    input: 
        # setup=expand('{tooling}/setup_status/rockfish.status', tooling=config['tooling_dir']),
        pod5=config['pod5_dir'] + "/{experiment}_{sr}kHz",
        bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: 
        main=directory("rockfish/{experiment}_{sr}kHz_{acc}_v{ver}"),
        pred="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}/predictions.tsv",
    threads: 40
    priority: 1
    resources: gpu=1
    log: "log/rockfish_call/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    conda: config['tooling']['rockfish']['conda']
    params: 
        model=config['tooling']['rockfish']['model'],
        ref=getRef
    shell: ntsh('''
        p5=$(realpath {input.pod5});
        mkdir -p {output.main}; 
        cd {output.main}; 
        {TIME} {ROCKFISH} inference \
            -i $p5 \
            --bam_path ../../{input.bam} \
            --model_path {params.model} \
            -r -t {threads} \
            -b 8192 -d 0;
        ''')

rule rockfish_map:
    input:  "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv"
    threads: 40
    priority: 1
    resources: gpu=1
    log: "log/rockfish_map/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    conda: config['tooling']['rockfish']['conda']
    params:
        ref=getRef
    shell: sh("python workflow/scripts/rockfish_extract_ref_pos.py --workers {threads} {input} {params.ref} > {output}")

rule rockfish_intersect:
    input: 
        ref="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv",
        pre="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}/predictions.tsv"
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_remapped.tsv"
    threads: 20
    log: "log/rockfish_intersect/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    priority: 1
    shell: sh("python workflow/scripts/rockfish_intersect.py -i {input.pre} -r {input.ref} -o {output}")

rule rockfish_aggregate:
    input:  "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_remapped.tsv"
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated.tsv"
    threads: 20
    log: "log/rockfish_aggregate/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("python workflow/scripts/rockfish_aggregate.py -i {input} -o {output}")

rule rockfish_getfasta:
    input:  
        bed="rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: "rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_consolidated_reference.tsv"
    threads: 20
    log: "log/rockfish_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: ntsh('''
        tmp=$(mktemp /tmp/rockfish_bed_fasta.XXXX);
        {BEDTOOLS} getfasta -fi {input.fasta} -bed {input.bed} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > $tmp;
        {TIME} {BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -s -i $tmp | {BEDTOOLS} getfasta -s -fi {input.fasta} -bed - -tab | cut -f 2 | paste $tmp - > {output};
        rm $tmp;
    ''')

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
