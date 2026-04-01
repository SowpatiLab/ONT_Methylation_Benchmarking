 
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

rule setup_deepmod2:
    output: config['tooling_dir'] + "/" + config['toolConfig']['deepmod2']['install_dir'] + "/" + config['toolConfig']['deepmod2']['executable']
    params: 
        conda_env=config['toolConfig']['deepmod2']['conda'],
        install_dir= config['tooling_dir'] + "/" + config['toolConfig']['deepmod2']['install_dir']
    container: None
    shell: """
        source "$(conda info --base)/etc/profile.d/conda.sh"

        [ ! -d {params.install_dir} ] && mkdir {params.install_dir} -p;
        cd {params.install_dir};

        if [[ -z "$(find . -maxdepth 0 -type d -empty 2>/dev/null)" ]];
        then
              echo "deepmod2 git already exists";
        else
            git clone https://github.com/WGLab/DeepMod2.git . ;
            cat src/detect.py \
                | perl -pe 's/^(from ont_fast5_api.*)/#\$1/' > tmp.txt;
            cat tmp.txt > src/detect.py;
            rm  tmp.txt;
        fi

        conda env create -f ./environment.yml -n {params.conda_env} -y;
        conda activate {params.conda_env};
        
        cuda_version=$(nvidia-smi | grep -oP 'CUDA Version: \\K[\\d.]+' | sed 's/\\.//')
        pip install torch torchvision --index-url https://download.pytorch.org/whl/cu${{cuda_version}}
        cd -
    """

rule deepmod2_call:
    input:  
        setup=config['tooling_dir'] + "/" + config['toolConfig']['deepmod2']['install_dir'] + "/" + config['toolConfig']['deepmod2']['executable'],
        pod5=config['pod5dir'] + "/{experiment}_{sr}kHz",
        bam=config['output_dir'] + "/" + "bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.cleansed.bam",
    output: directory(config['output_dir'] + "/" + "tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.cleansed")
    threads: 40
    priority: 2
    resources: gpu=1
    container: None
    # conda: config['toolConfig']['deepmod2']['conda']
    log: "log/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.log"
    params: 
        ref=getRef,
        model=inferDeepMod2Model,
        exec=config['tooling_dir'] + "/" + config['toolConfig']['deepmod2']['install_dir'] + "/" + config['toolConfig']['deepmod2']['executable'],
        conda_env=config['toolConfig']['deepmod2']['conda']
    shell:  ntsh("""
        exec="conda run -n {params.conda_env} python {params.exec}";
        {TIME} ${{exec}} detect \
            --bam {input.bam} \
            --input {input.pod5} \
            --model {params.model} \
            --file_type pod5 \
            --threads {threads} \
            --ref {params.ref} \
            --seq_type dna \
            --batch_size 2048 \
            --output {output} --device=cuda""")

rule deepmod_rebed_add_ref:
    input:   
        bed=config['output_dir'] + "/" + "tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.cleansed",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output: config['output_dir'] + "/" + "tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.tsv"
    log:     "log/deepmod2_addref/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.log"
    threads: 20
    priority: 0
    params: 
        ref=getRef
    conda: str(workflow.basedir) + "/" + config['default_conda_env']
    shell:   ntsh('''
        bed=$(mktemp /tmp/deepmod2_metadata.XXXX);
        sed '1d' {input.bed}/output.per_site | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$6,$5,$4,$7,$8,$9}}' > $bed;
        
        {TIME} {BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i $bed -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $bed - > {output} 
    ''')

rule colsolidate_deepmod2:
    input:  config['output_dir'] + "/" + "tool_out/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.tsv"
    output: config['output_dir'] + "/" + "meta/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.std.bed"
    threads: 20
    priority: 0
    conda: str(workflow.basedir) + "/" + config['default_conda_env']
    log: "log/colsolidate_deepmod2_{deepmodel}/{experiment}_{sr}kHz_{acc}_v{ver}.log"
    params: script_dir=config['scripts_common']
    shell: sh("python {params.script_dir}/deepmod2_consolidate.py {input} {output} {wildcards.deepmodel}")
