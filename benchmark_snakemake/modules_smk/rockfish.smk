rule setup_rockfish:
    output: config['output_dir'] + "/" + "{tooling}/setup_status/rockfish.status"
    params: 
        conda_env=config['toolConfig']['rockfish']['conda'],
        model=subpath(config['toolConfig']['rockfish']['model'], parent=True)
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

rule download_rockfish_model:
    output: ROCKFISH_MODEL
    conda: config['toolConfig']['rockfish']['conda']
    params:
        save_dir=subpath(output[0], parent=True),
        script_dir=config['scripts_common']
    shell: """
        [ ! -d {params.save_dir} ] && mkdir -p {params.save_dir} ;
        micromamba run -n rockfish python {params.script_dir}/rockfish_model_dl.py download -m 5kHz -s {params.save_dir}
    """
    ## the original script was edited since cookies were causing issues
    # rockfish download -m 5kHz -s {params.save_dir} ;

rule rockfish_call:
    input: 
        # setup=expand('{tooling}/setup_status/rockfish.status', tooling=config['tooling_dir']),
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz",
        bam=config['output_dir'] + "/" + "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam",
        model=ROCKFISH_MODEL
    output: 
        main=directory(config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}"),
        pred=config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}/predictions.tsv",
    threads: 40
    priority: 1
    resources: gpu=1
    log: "log/rockfish_call/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    conda: config['toolConfig']['rockfish']['conda']
    params: 
        # model=ROCKFISH_MODEL,
        ref=getRef,
        call_flags=config['toolConfig']['rockfish']['call_flags'],
    shell:  ntsh("""
        p5=$(realpath {input.pod5});
        bam_dir=$(realpath {input.bam});
        model=$(realpath {input.model});
        
        mkdir -p {output.main}; 
        cd {output.main}; 
        {TIME} rockfish inference \
            -i $p5 \
            --bam_path $bam_dir \
            --model_path $model \
            -t {threads} {params.call_flags};
    """)

rule rockfish_map:
    input:  config['output_dir'] + "/" + "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam",
    output: config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv"
    threads: 40
    priority: 1
    resources: gpu=1
    log: "log/rockfish_map/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    conda: config['toolConfig']['rockfish']['conda']
    params:
        ref=getRef,
        script_dir=config['scripts_common']
    shell: ntsh("""
        if  [ -f /.dockerenv ];
        then
            eval "$(micromamba shell hook --shell bash)";
            micromamba activate rockfish;
        elif  [ -d "/.singularity.d" ];
        then
            eval "$(micromamba shell hook --shell bash)";
            micromamba activate rockfish;
        fi;
        {TIME} python {params.script_dir}/rockfish_extract_ref_pos.py \
            --workers {threads} {input} {params.ref} > {output}""")

rule rockfish_intersect:
    input: 
        ref=config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}_bam_ref_map.tsv",
        pre=config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}/predictions.tsv"
    output: config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.remapped.tsv"
    threads: 20
    log: "log/rockfish_intersect/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    priority: 1
    params: 
        script_dir=config['scripts_common']
    shell: sh("python {params.script_dir}/rockfish_intersect.py -i {input.pre} -r {input.ref} -o {output}")

rule rockfish_aggregate:
    input:  config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.remapped.tsv"
    output: config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.aggregated.tsv"
    threads: 20
    log: "log/rockfish_aggregate/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        script_dir=config['scripts_common']
    shell: sh("python {params.script_dir}/rockfish_aggregate.py -i {input} -o {output}")

rule rockfish_getfasta:
    input:  
        bed=config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.aggregated.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.aggregated.rebed.ref.tsv"
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
    input:  config['output_dir'] + "/" + "tool_out/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.aggregated.rebed.ref.tsv"
    output: config['output_dir'] + "/" + "meta/rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.rockfish.aggregated.rebed.ref.std.bed"
    threads: 20
    conda: str(workflow.basedir) + "/" + config['default_conda_env']
    log: "log/consolidate_rockfish/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        script_dir=config['scripts_common']
    shell: sh("python {params.script_dir}/rockfish_consolidate.py {input} {output}")
