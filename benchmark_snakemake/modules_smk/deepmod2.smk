 
def inferDeepMod2Model(wildcards):
    models = {
        '5kHz_transformer_v5.2r2': "transformer_r10.4.1_5khz_v5.0",
        '5kHz_transformer_v5.2r1': "transformer_r10.4.1_5khz_v5.0",
        '5kHz_transformer_v5r3': "transformer_r10.4.1_5khz_v5.0",
        '5kHz_transformer_v4r2': "transformer_r10.4.1_5khz_v4.3",
        '4kHz_transformer_v4r2': "transformer_r10.4.1_4khz_v4.1",
        '5kHz_BiLSTM_v5.2r2'     : "bilstm_r10.4.1_5khz_v5.0",
        '5kHz_BiLSTM_v5.2r1'     : "bilstm_r10.4.1_5khz_v5.0",
        '5kHz_BiLSTM_v5r3'     : "bilstm_r10.4.1_5khz_v5.0",
        '5kHz_BiLSTM_v4r2'     : "bilstm_r10.4.1_5khz_v4.3",
        '4kHz_BiLSTM_v4r2'     : "bilstm_r10.4.1_4khz_v4.1",
    }
    sr, arch, ver = wildcards.sr, wildcards.deepmodel, wildcards.ver
    return models[f'{sr}kHz_{arch}_v{ver}']

# rule setup_deepmod2:
#     output: "{tooling}/setup_status/deepmod2.status"
#     params: 
#         conda_env=config['toolConfig']['deepmod2']['conda']
#     shell: """
#         source "$(conda info --base)/etc/profile.d/conda.sh"
#         [ ! -d {wildcards.tooling} ] && mkdir {wildcards.tooling} -p
#         cd {wildcards.tooling}

#         git clone https://github.com/WGLab/DeepMod2.git deepmod2
#         conda env create -f deepmod2/environment.yml -n {params.conda_env} -y
#         conda activate {params.conda_env}
        
#         cuda_version=$(nvidia-smi | grep -oP 'CUDA Version: \\K[\\d.]+' | sed 's/\\.//')
#         echo $cuda_version
#         pip install torch torchvision --index-url https://download.pytorch.org/whl/cu${{cuda_version}}
#         cd ../
#         echo deepmod2 dir: $(realpath deepmod2) >> {output}
#     """

rule deepmod2_call:
    input:  
        # setup=expand("{tooling}/setup_status/deepmod2.status", tooling=config['tooling_dir']),
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz",
        bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam",
    output: directory("tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.cleansed")
    threads: 40
    priority: 2
    resources: gpu=1
    conda: config['toolConfig']['deepmod2']['conda']
    log: "log/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.log"
    params: 
        ref=getRef,
        model=inferDeepMod2Model
    shell: sh("python {DEEPMOD2} detect \
                --bam {input.bam} \
                --input {input.pod5} \
                --model {params.model} \
                --file_type pod5 \
                --threads {threads} \
                --ref {params.ref} \
                --seq_type dna \
                --batch_size 2048 \
                --output {output} --device=cuda")

rule deepmod_rebed_add_ref:
    input:   
        bed="tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.cleansed",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output:  "tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.tsv"
    log:     "log/deepmod2_addref/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.log"
    threads: 20
    priority: 0
    params: 
        ref=getRef
    conda:  f"../{config['std_conda']}"
    shell:   ntsh('''
        bed=$(mktemp /tmp/deepmod2_metadata.XXXX);
        sed '1d' {input.bed}/output.per_site | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$6,$5,$4,$7,$8,$9}}' > $bed;
        
        {TIME} {BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i $bed -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $bed - > {output} 
    ''')

rule colsolidate_deepmod2:
    input:  "tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.tsv"
    output: "meta/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.std.bed"
    threads: 20
    priority: 0
    conda:  f"../{config['std_conda']}"
    log: "log/colsolidate_deepmod2_{deepmodel}/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    shell: sh("python workflow/scripts_common/deepmod2_consolidate.py {input} {output} {wildcards.deepmodel}")
