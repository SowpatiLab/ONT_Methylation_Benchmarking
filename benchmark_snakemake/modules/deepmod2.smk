
def inferDeepMod2Model(wildcards):
    models = {
        '5kHz_transformer_v5r3': "transformer_r10.4.1_5khz_v5.0",
        '5kHz_transformer_v4r2': "transformer_r10.4.1_5khz_v4.3",
        '4kHz_transformer_v4r2': "transformer_r10.4.1_4khz_v4.1",
        '5kHz_BiLSTM_v5r3'     : "bilstm_r10.4.1_5khz_v5.0",
        '5kHz_BiLSTM_v4r2'     : "bilstm_r10.4.1_5khz_v4.3",
        '4kHz_BiLSTM_v4r2'     : "bilstm_r10.4.1_4khz_v4.1",
    }
    sr, arch, ver = wildcards.sr, wildcards.deepmodel, wildcards.ver
    return models[f'{sr}kHz_{arch}_v{ver}']

rule setup_deepmod2:
    output: "{tooling}/setup_status/deepmod2.status"
    params: 
        conda_env=config['tooling']['deepmod2']['conda']
    shell: """
        source "$(conda info --base)/etc/profile.d/conda.sh"
        [ ! -d {wildcards.tooling} ] && mkdir {wildcards.tooling} -p
        cd {wildcards.tooling}

        git clone https://github.com/WGLab/DeepMod2.git deepmod2
        conda env create -f deepmod2/environment.yml -n {params.conda_env} -y
        conda activate {params.conda_env}
        
        cuda_version=$(nvidia-smi | grep -oP 'CUDA Version: \\K[\\d.]+' | sed 's/\\.//')
        echo $cuda_version
        pip install torch torchvision --index-url https://download.pytorch.org/whl/cu${{cuda_version}}
        cd ../
        echo deepmod2 dir: $(realpath deepmod2) >> {output}
    """

rule deepmod2_call:
    input:  
        # setup=expand("{tooling}/setup_status/deepmod2.status", tooling=config['tooling_dir']),
        pod5=config['pod5_dir'] + "/{experiment}_{sr}kHz",
        bam="bam/sorted_move_cleansed/{experiment}_{sr}kHz_{acc}_v{ver}.bam",
    output: directory("deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}")
    threads: 40
    priority: 2
    resources: gpu=1
    conda: config['tooling']['deepmod2']['conda']
    log: "log/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.log"
    params: 
        ref=getRef,
        model=inferDeepMod2Model,
        execut=lambda x: f"{config['tooling_dir']}/{config['tooling']['deepmod2']['install']}/{config['tooling']['deepmod2']['execut']}"
    shell: sh("python {params.execut} detect \
                --bam {input.bam} \
                --input {input.pod5} \
                --model {params.model} \
                --file_type pod5 \
                --threads {threads} \
                --ref {params.ref} \
                --seq_type dna \
                --batch_size 2048 \
                --output {output} --device=cuda")

rule consolidate_ref:
    input:   
        bed="deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}",
        fasta=lambda wildcards: getRef(wildcards),
        genome=lambda wildcards: re.sub(r'.fa(sta|)', '.genome', getRef(wildcards))
    output:  "deepmod2/ref/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.tsv"
    log:     "log/deepmod2_addref/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.tsv"
    threads: 20
    priority: 0
    params: 
        ref=getRef
    shell:   ntsh('''
        bed=$(mktemp /tmp/deepmod2_metadata.XXXX);
        sed '1d' {input.bed}/output.per_site | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$6,$5,$4,$7,$8,$9}}' > $bed;
        
        {TIME} {BEDTOOLS} slop -g {input.genome} -l 5 -r 11 -i $bed -s | {BEDTOOLS} getfasta -fi {input.fasta} -bed - -tab -s | cut -f 2 | paste $bed - > {output} 
    ''')

rule colsolidate_deepmod2:
    input:  "deepmod2/ref/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.tsv"
    output: "meta/deepmod2/{experiment}_{sr}kHz_{acc}_v{ver}_{deepmodel}.tsv"
    threads: 20
    priority: 0
    run:
        import polars as pl
        from pathlib import Path
        import re

        t = pl.read_csv(f"{input[0]}", separator='\t', has_header=True, new_columns=['chrom', 'p1', 'p2', 'coverage', 'is_cpg', 'strand', 'M', 'UM', 'per', 'full_context'])
        select  = ['chrom', 'p1', 'p2', 'mod', 'coverage', 'strand', 'M', 'UM', 'per', 'flowcell', 'tool', 'model', 'sample', 'acc', 'sample_rate', 'species', 'full_context']
        species = wildcards.experiment.split('_')[0]
        
        (
            t.with_columns([
                pl.lit('5mC').alias('mod'),
                pl.lit('r10.4.1').alias('flowcell'),
                pl.lit(wildcards.experiment).alias('sample'),
                pl.lit(f'deepmod2_{wildcards.deepmodel}').alias('tool'),
                pl.lit(wildcards.acc).alias('acc'),
                pl.lit(wildcards.deepmodel).alias('model'),
                pl.lit(f"{wildcards.sr}kHz").alias('sample_rate'),
                pl.lit(species).alias('species'),
                (pl.col('per')*100).alias('per')
            ])
            .select(select)
        ).write_csv(output[0], separator='\t', include_header=True)
