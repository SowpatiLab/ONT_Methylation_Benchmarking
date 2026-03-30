
# rule setup_f5c:
#     output: config['output_dir'] + "/" + "{tooling}/setup_status/f5c_{version}.status"
#     params: version="v1.5"
#     shell: """
#         [ ! -d {wildcards.tooling} ] && mkdir {wildcards.tooling} -p
#         cd {wildcards.tooling}
        
#         wget "https://github.com/hasindu2008/f5c/releases/download/{params.version}/f5c-{params.version}-binaries.tar.gz" && \
#             tar xvf f5c-{params.version}-binaries.tar.gz
        
#         cd ../
#         echo "f5c dir: $(realpath {wildcards.tooling}/f5c-{params.version})" > {output}
#     """

rule pod5_to_blow5:
    input:  config['pod5dir'] + "/{experiment}_{sr}kHz"
    output: config['output_dir'] + "/" + "intermediary/blow5/{experiment}_{sr}kHz.blow5"
    threads: 20
    log: "log/pod5_to_blow5/{experiment}_{sr}kHz.log"
    shell: ntsh("{TIME} blue-crab p2s {input} -t {threads} -o {output}")

rule f5c_idx_fastq:
    input:  
        # setup=expand('{tooling}/setup_status/f5c_{version}.status', tooling=config['tooling_dir'], version=config['toolConfig']['f5c']['version']),
        blow5=config['output_dir'] + "/" + "intermediary/blow5/{experiment}_{sr}kHz.blow5",
        fastq=config['output_dir'] + "/" + "intermediary/fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq"
    output: 
        fastq=multiext(config['output_dir'] + "/" + "intermediary/fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq", '.index', '.index.fai', '.index.gzi')
    threads: 20
    log: "log/f5c_idx/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params:
        exec=config['toolConfig']['f5c']['executable']
    shell:  ntsh("""
        if  [ -f /.dockerenv ];
        then
            exec={params.exec};
        elif [ -d "/.singularity.d" ];
        then
            exec={params.exec};
        else
            exec={F5C};
        fi;

        {TIME} ${{exec}} index --slow5 {input.blow5} {input.fastq};""")

rule f5c_call:
    input: 
        # setup=expand('{tooling}/setup_status/f5c_{version}.status', tooling=config['tooling_dir'], version=config['toolConfig']['f5c']['version']),
        blow5=config['output_dir'] + "/" + "intermediary/blow5/{experiment}_{sr}kHz.blow5",
        bam=config['output_dir'] + "/" + "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam",
        fastq=config['output_dir'] + "/" + "intermediary/fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq",
        fastq_idx=config['output_dir'] + "/" + "intermediary/fastq/{experiment}_{sr}kHz_{acc}_v{ver}.fastq.index"
    output: config['output_dir'] + "/" + "tool_out/f5c/readwise/{experiment}_{sr}kHz_{acc}_v{ver}.f5c.tsv"
    threads: 20
    priority: 1
    resources: gpu=1
    log: "log/f5c_call/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: 
        ref=getRef,
        call_flags=config['toolConfig']['f5c']['call_flags'],
        exec=config['toolConfig']['f5c']['executable']
    shell:  ntsh("""
        if  [ -f /.dockerenv ];
        then
            exec={params.exec};
        elif [ -d "/.singularity.d" ];
        then
            exec={params.exec};
        else
            exec={F5C};
        fi;

        {TIME} ${{exec}} call-methylation \
            {params.call_flags} \
            -t {threads} \
            --slow5 {input.blow5} \
            -b {input.bam} \
            -g {params.ref} \
            -r {input.fastq} > {output}""")

rule f5c_restrand:
    input:  config['output_dir'] + "/" + "tool_out/f5c/readwise/{experiment}_5kHz_{acc}_v{ver}.f5c.tsv"
    output: config['output_dir'] + "/" + "tool_out/f5c/readwise_stranded/{experiment}_5kHz_{acc}_v{ver}.f5c.stranded.tsv"
    log:    "log/f5c_restrand/{experiment}_5kHz_{acc}_v{ver}.tsv"
    threads: 20
    params: 
        ref=getRef,
        script_dir=config['scripts_common']
    shell: sh("python {params.script_dir}/f5c_restrand.py -i {input} -r {params.ref} -o {output}")

rule f5c_aggregate:
    input:  
        # setup=expand('{tooling}/setup_status/f5c_{version}.status', tooling=config['tooling_dir'], version=config['toolConfig']['f5c']['version']),
        out=config['output_dir'] + "/" + "tool_out/f5c/readwise{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.tsv"
    output: config['output_dir'] + "/" + "tool_out/f5c/aggregated{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.aggregated.tsv"
    log: "log/f5c_aggregate{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}{stranded_ext}.log"
    threads: 20
    params:
        exec=config['toolConfig']['f5c']['executable']
    shell: ntsh("""
        if  [ -f /.dockerenv ];
        then
            exec={params.exec};
        elif [ -d "/.singularity.d" ];
        then
            exec={params.exec};
        else
            exec={F5C};
        fi;

        {TIME} ${{exec}} meth-freq -i {input.out} -s > {output};
    """)

rule f5c_fasta:
    input:  
        bed=config['output_dir'] + "/" + "tool_out/f5c/aggregated{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.aggregated.tsv",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: config['output_dir'] + "/" + "tool_out/f5c/aggregated_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.aggregated.rebed.ref.tsv"
    log:    "log/f5c_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.log"
    threads: 20
    run: 
        if wildcards.strandstate=="_stranded":
            shell(ntsh('''
                tmp=$(mktemp /tmp/f5c_bed_fasta.XXXX);
                
                {BEDTOOLS} getfasta -fi {input.fasta} -bed <(tail -n +2 {input.bed}) -bedOut | \
                    awk 'BEGIN{{OFS="\\t"}} {{ if(toupper($9)=="C") $9="+"; else if(toupper($9)=="G") $9="-"; print $1,$2,$3,$5,".",$9,$6,$5-$6,$7 }}' > $tmp;

                {TIME} {BEDTOOLS} slop -l 5 -r 11 -s -g {input.genome} -i $tmp \
                    | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $tmp - > {output} {STDERR};
                rm $tmp
            '''))
        else:
            shell(ntsh('''
                tmp=$(mktemp /tmp/f5c_bed_fasta.XXXX);

                sed '1d' {input.bed} | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$2+1,$5,".","+",$6,$5-$6,$7}}' > $tmp;
                {TIME} {BEDTOOLS} slop -l 5 -r 11 -s -g {input.genome} -i $tmp \
                    | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $tmp - > {output} {STDERR};
                rm $tmp
            '''))
            
rule consolidate_f5c:
    input:  config['output_dir'] + "/" + "tool_out/f5c/aggregated_fasta{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.aggregated.rebed.ref.tsv"
    output: config['output_dir'] + "/" + "meta/f5c{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}.f5c{stranded_ext}.aggregated.rebed.ref.std.bed"
    threads: 20
    log: "log/consolidate_f5c{strandstate}/{experiment}_{sr}kHz_{acc}_v{ver}{stranded_ext}.log"
    conda: str(workflow.basedir) + "/" + config['default_conda_env']
    params: script_dir=config['scripts_common']
    shell: sh("python {params.script_dir}/f5c_consolidate.py {input} {output}")