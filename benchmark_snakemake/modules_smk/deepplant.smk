rule bam_fn_reorder:
    input:  config['output_dir'] + "/" + "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam"
    output: config['output_dir'] + "/" + "bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    log:    "log/sorted_move_by_fntag/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    conda:  config['default_conda_env']
    threads: 20
    resources:
        queue='hm-q',
        mem='250gb'
    shell:  sh("{SAMTOOLS} sort -t fn {input} -@ {threads} -O BAM -o {output}")

rule download_deepplant_model:
    output: directory(DEEPPLANT_MODEL)
    params:
        tooling_dir=config['tooling_dir'],
        install_dir=config['toolConfig']['deepplant']['install_dir'],
        save_dir=subpath(output[0], parent=True)
    shell: """
        git clone https://github.com/xiaochuanle/DeepPlant.git {params.tooling_dir}/{params.install_dir};
        [ ! -d {params.save_dir} ] && mkdir -p {params.save_dir};
        cp -r {params.tooling_dir}/{params.install_dir}/model/* {params.save_dir};
    """


rule deepplant_call:
    input:  
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz",
        bam=config['output_dir'] + "/" + "bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
        model=DEEPPLANT_MODEL
    output: 
        cpg=config['output_dir'] + "/" + "tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_cpg.tsv",
        chg=config['output_dir'] + "/" + "tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_chg.tsv",
        chh=config['output_dir'] + "/" + "tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_chh.tsv"
    threads: 32
    log: "log/deeplant/{experiment}_{sr}kHz_{acc}_v{ver}_v{ver}.log"
    conda: config['toolConfig']['deepplant']['conda']
    priority: 3
    params: 
        ref=getRef,
        call_flags=config['toolConfig']['deepplant']['call_flags'],
        parent_dir=subpath(output[0], parent=True),
        exp_name=subpath(output[0], strip_suffix=".deepplant_cpg.tsv"),
        # exec=config['toolConfig']['deepplant']['executable']
    shell: ntsh("""
        [ ! -d {params.exp_name} ] && mkdir -p {params.exp_name};

        {TIME} {DEEPPLANT} extract_and_call_mods \
            {input.pod5} {input.bam} {params.ref} DNA {params.exp_name} {input.model} \
            {params.call_flags};
        
        mv {params.exp_name}/cpg_result.txt {params.exp_name}.deepplant_cpg.tsv;
        mv {params.exp_name}/chg_result.txt {params.exp_name}.deepplant_chg.tsv;
        mv {params.exp_name}/chh_result.txt {params.exp_name}.deepplant_chh.tsv;

        rm -r {params.exp_name};
    """)

rule deepplant_aggregate:
    input:  config['output_dir'] + "/" + "tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.tsv"
    output: config['output_dir'] + "/" + "tool_out/deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.tsv"
    threads: 20
    conda:  config['default_conda_env']
    log:    "log/deepplant_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: 
        aggregation_flags=config['toolConfig']['deepbam']['aggregation_flags'],
        script_dir=Path(workflow.basedir) / "scripts_common"
    shell:  sh("python {params.script_dir}/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} {params.aggregation_flags}")

rule deepplant_rebed:
    input:  config['output_dir'] + "/" + "tool_out/deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.tsv"
    output: config['output_dir'] + "/" + "tool_out/deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.tsv"
    threads: 20
    log:    "log/deepplant_rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: ref=getRef
    shell: sh('''
        {BEDTOOLS} getfasta -fi {params.ref} -bed {input} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output};
    ''')

rule deepplant_addfasta:
    input: 
        bed=config['output_dir'] + "/" + "tool_out/deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: config['output_dir'] + "/" + "tool_out/deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.ref.tsv"
    threads: 20
    log:    "log/deepplant_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell: sh("{BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i {input.bed} -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste {input.bed} - > {output}")

rule consolidate_deepplant:
    input:  config['output_dir'] + "/" + "tool_out/deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.ref.tsv"
    output: config['output_dir'] + "/" + "meta/deepplant/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.ref.std.bed"
    threads: 20
    conda:  config['default_conda_env']
    log: "log/deepplant/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: script_dir=Path(workflow.basedir) / "scripts_common"
    shell: sh("python {params.script_dir}/deeptools_consolidate.py {input} {output}")
    # shell: "python {workflow.basedir}/scripts_common/deeptools_consolidate.py {input} {output}"

