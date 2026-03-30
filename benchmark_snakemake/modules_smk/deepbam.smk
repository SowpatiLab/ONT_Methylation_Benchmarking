
rule download_deepbam_model:
    output: DEEPBAM_MODEL
    params:
        tooling_dir=config['tooling_dir'],
        install_dir=config['toolConfig']['deepbam']['install_dir'],
        save_dir=subpath(output[0], parent=True)
    shell: """
        if [ -d {params.tooling_dir}/{params.install_dir} ];
        then
            echo {params.tooling_dir}/{params.install_dir} exists;
        elif [ -d {params.tooling_dir}/{params.install_dir}/traced_script_module/ ];
        then
            echo {params.tooling_dir}/{params.install_dir}/traced_script_module/ exists;
        else 
            git clone https://github.com/xiaochuanle/DeepBAM.git {params.tooling_dir}/{params.install_dir};
        fi

        [ ! -d {params.save_dir} ] && mkdir -p {params.save_dir};
        ls {params.tooling_dir}/{params.install_dir}/traced_script_module/;
        cp -r {params.tooling_dir}/{params.install_dir}/traced_script_module/* {params.save_dir}
    """

rule deepbam_call:
    input:  
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz/",
        bam=config['output_dir'] + "/" + "bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
        model=DEEPBAM_MODEL
    output: config['output_dir'] + "/" + "tool_out/deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.tsv"
    threads: 32
    priority: 0
    resources: gpu=1
    conda: config['toolConfig']['deepbam']['conda']
    log: "log/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        ref=getRef,
        # model=DEEPBAM_MODEL,
        call_flags=config['toolConfig']['deepbam']['call_flags']
    shell: ntsh('''
        DeepBAM extract_and_call_mods \
        {input.pod5} {input.bam} {params.ref} \
        DNA {output} {input.model} {params.call_flags}
    ''')

rule deepbam_aggregate:
    input:  config['output_dir'] + "/" + "tool_out/deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.tsv"
    output: config['output_dir'] + "/" + "tool_out/deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.tsv"
    threads: 20
    log:    "log/deepbam_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        aggregation_flags=config['toolConfig']['deepbam']['aggregation_flags'],
        script_dir=Path(f"{workflow.basedir}/scripts_common") 
    shell:  sh("python {params.script_dir}/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} {params.aggregation_flags}")

rule deepbam_rebed:
    input:  
        bed=config['output_dir'] + "/" + "tool_out/deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.tsv",
        fasta=lambda wildcards: getRef(wildcards)
    output: config['output_dir'] + "/" + "tool_out/deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.tsv"
    threads: 20
    log:    "log/deepbam_rebed/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: ref=getRef
    shell: sh('''
        {BEDTOOLS} getfasta -fi {input.fasta} -bed {input.bed} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output};
    ''')

rule deepbam_addfasta:
    input:  
        bed=config['output_dir'] + "/" + "tool_out/deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: config['output_dir'] + "/" + "tool_out/deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.ref.tsv"
    threads: 20
    log:    "log/deepbam_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("{BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i {input.bed} -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste {input.bed} - > {output}")

rule consolidate_deepbam:
    input:  config['output_dir'] + "/" + "tool_out/deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.ref.tsv"
    output: config['output_dir'] + "/" + "meta/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.ref.std.bed"
    threads: 20
    conda: str(workflow.basedir) + "/" + config['default_conda_env']
    log: "log/consolidate_deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: script_dir=Path(f"{workflow.basedir}/scripts_common") 
    shell: sh("python {params.script_dir}/deeptools_consolidate.py {input} {output}")
