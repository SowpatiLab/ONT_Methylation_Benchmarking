rule bam_fn_reorder:
    input:  "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam"
    output: "bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    log:    "log/sorted_move_by_fntag/{experiment}_{sr}kHz_{acc}_v{ver}.bam"
    conda:  f"../{config['std_conda']}"
    threads: 20
    resources:
        queue='hm-q',
        mem='250gb'
    shell:  sh("{SAMTOOLS} sort -t fn {input} -@ {threads} -O BAM -o {output}")

rule deepplant_call:
    input:  
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz",
        bam="bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: 
        cpg="tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_cpg.tsv",
        chg="tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_chg.tsv",
        chh="tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_chh.tsv"
    threads: 32
    log: "log/deeplant/{experiment}_{sr}kHz_{acc}_v{ver}_v{ver}.log"
    conda: config['toolConfig']['deepplant']['conda']
    priority: 3
    params: 
        ref=getRef,
        model_dir=config['toolConfig']['deepplant']['model'],
        call_flags=config['toolConfig']['deepplant']['call_flags'],
        parent_dir=subpath(output[0], parent=True),
        exp_name=subpath(output[0], strip_suffix=".deepplant_cpg.tsv"),
    shell: ntsh('''
        [ ! -d {params.exp_name} ] && mkdir -p {params.exp_name};

        {TIME} {DEEPPLANT} extract_and_call_mods \
            {input.pod5} {input.bam} {params.ref} DNA {params.exp_name} {params.model_dir} \
            {params.call_flags};
        
        mv {params.exp_name}/cpg_result.txt {params.exp_name}.deepplant_cpg.tsv;
        mv {params.exp_name}/chg_result.txt {params.exp_name}.deepplant_chg.tsv;
        mv {params.exp_name}/chh_result.txt {params.exp_name}.deepplant_chh.tsv;

        rm -r {params.exp_name};
    ''')

rule deepplant_aggregate:
    input:  "tool_out/deepplant/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.tsv"
    output: "tool_out/deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.tsv"
    threads: 20
    conda:  f"../{config['std_conda']}"
    log:    "log/deepplant_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: aggregation_flags=config['toolConfig']['deepbam']['aggregation_flags']
    shell:  sh("python workflow/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} {params.aggregation_flags}")

rule deepplant_rebed:
    input:  "tool_out/deepplant/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.tsv"
    output: "tool_out/deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.tsv"
    threads: 20
    log:    "log/deepplant_rebed/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    params: ref=getRef
    shell: sh('''
        {BEDTOOLS} getfasta -fi {params.ref} -bed {input} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output};
    ''')

rule deepplant_addfasta:
    input: 
        bed="tool_out/deepplant/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: "tool_out/deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.ref.tsv"
    threads: 20
    log:    "log/deepplant_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell: sh("{BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i {input.bed} -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste {input.bed} - > {output}")

rule consolidate_deepplant:
    input:  "tool_out/deepplant/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.ref.tsv"
    output: "meta/deepplant/{experiment}_{sr}kHz_{acc}_v{ver}.deepplant_{context}.aggregated.rebed.ref.std.bed"
    threads: 20
    conda:  f"../{config['std_conda']}"
    log: "log/deepplant/{experiment}_{sr}kHz_{acc}_v{ver}_{context}.log"
    shell: sh("python workflow/scripts_common/deeptools_consolidate.py {input} {output}")
    # shell: "python {workflow.basedir}/scripts_common/deeptools_consolidate.py {input} {output}"

