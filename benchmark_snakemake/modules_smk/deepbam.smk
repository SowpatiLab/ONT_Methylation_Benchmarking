rule deepbam_call:
    input:  
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz/",
        bam="bam/sorted_move_fnord/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: "tool_out/deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.tsv"
    threads: 32
    priority: 0
    resources: gpu=1
    conda: config['toolConfig']['deepbam']['conda']
    container: "/data1/ccmb/reuben/benchmarking/bacterial_nf/apptainer_final/ontMethylationBenchmarking.sif"
    log: "log/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        ref=getRef,
        model=config['toolConfig']['deepbam']['model'],
        call_flags=config['toolConfig']['deepbam']['call_flags']
    shell: ntsh('''
        DeepBAM extract_and_call_mods \
        {input.pod5} {input.bam} {params.ref} \
        DNA {output} {params.model} {params.call_flags}
    ''')

rule deepbam_aggregate:
    input:  "tool_out/deepbam/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.tsv"
    output: "tool_out/deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.tsv"
    threads: 20
    log:    "log/deepbam_aggregation/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: aggregation_flags=config['toolConfig']['deepbam']['aggregation_flags']
    shell:  sh("python workflow/deepbam_aggregate.py \
                --file_path  {input} \
                --aggregation_output {output} {params.aggregation_flags}")

rule deepbam_rebed:
    input:  
        bed="tool_out/deepbam/aggregated/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.tsv",
        fasta=lambda wildcards: getRef(wildcards)
    output: "tool_out/deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.tsv"
    threads: 20
    log:    "log/deepbam_rebed/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: ref=getRef
    shell: sh('''
        {BEDTOOLS} getfasta -fi {input.fasta} -bed {input.bed} -bedOut | \
            awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($7)=="C") $7="+"; else if(toupper($7)=="G") $7="-"; print $1,$2,$3,$5+$6,".",$7,$5,$6,$4 }}' > {output};
    ''')

# //todo: combine addref and rebed 
rule deepbam_addfasta:
    input:  
        bed="tool_out/deepbam/rebed/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: "tool_out/deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.ref.tsv"
    threads: 20
    log:    "log/deepbam_addfasta/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("{BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i {input.bed} -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste {input.bed} - > {output}")

rule consolidate_deepbam:
    input:  "tool_out/deepbam/agg_fasta/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.ref.tsv"
    output: "meta/deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.deepbam.aggregated.rebed.ref.std.bed"
    threads: 20
    conda:  f"../{config['std_conda']}"
    log: "log/consolidate_deepbam/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("python workflow/scripts_common/deeptools_consolidate.py {input} {output}")
