include {
    select_base_model ;
    select_mod_model ;
    fetch_ref ;
} from '../lib/utils.groovy'


process install_rerio {
    label 'download'
    storeDir "${params.tooling_dir}"

    output:
        path "rerio"
        // path "res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1",            emit: _4kHz_sup_v4r1
        // path "res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2",     emit: _4kHz_sup_v4r1_5mC
        // path "res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2",     emit: _4kHz_sup_v4r1_6mA
        // path "res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1", emit: _5kHz_sup_v4r1_4mC

    shell:
    """
        git clone https://github.com/nanoporetech/rerio.git rerio
    """
}

process download_rerio_models {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.model_dir}"

    input:  tuple val(base_model), val(key),  val(model_name), val(dorado), val(rerio)
    output: tuple val(base_model), val(key), path(model_name)

    shell:
    """
        python ${rerio}/download_model.py --dorado ${rerio}/dorado_models/${model_name}_url
        cp -r ${rerio}/dorado_models/${model_name} ./
    """

}

process install_dorado {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.install_dir}"

    input:
        val dorado_version_map

    output:
        path "${params.toolConfig.dorado.dorado_version_aliases[dorado_version_map]}"

    script:
    """
        [ ! -f ${dorado_version_map}-linux-x64.tar.gz ] && wget https://cdn.oxfordnanoportal.com/software/analysis/${dorado_version_map}-linux-x64.tar.gz \
            -O ${dorado_version_map}-linux-x64.tar.gz
        tar -xvzf ${dorado_version_map}-linux-x64.tar.gz 
        mv "${dorado_version_map}"-linux-x64 "${params.toolConfig.dorado.dorado_version_aliases[dorado_version_map]}"
    """
}

process download_dorado_base_model {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.model_dir}"

    input:
        tuple val(model_name), val(dorado)

    output:
        tuple val(model_name), val(dorado), path(model_name)

    shell:
    """
        ${dorado}/bin/dorado download --model ${model_name};
    """
}

process download_dorado_mod_models {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.model_dir}"

    input:
        tuple val(base_model), val(key), val(model_name), val(dorado)

    output:
        tuple val(base_model), val(key), path(model_name)

    shell:
    """
        ${dorado}/bin/dorado download --model ${model_name};
    """
}

process dorado_move {
    storeDir "bam/sorted_move"
    
    label 'gpu'
    label 'std_conda'
    // conda 'bioconda::samtools==1.21'

    input:
        tuple val(experiment), val(key), val(dorado), path(base_model)

    output:
    path "${experiment}_${key}.bam"
    path "${experiment}_${key}.bam.csi"

    script:
    reference = fetch_ref(experiment)
    """
        ${dorado}/bin/dorado basecaller \
            ${base_model} ${params.pod5dir}/${experiment}_5kHz \
            ${params.toolConfig.dorado.run_flags} --reference ${reference} \
            --emit-moves | samtools sort -O BAM -@ ${task.cpus} -o "${experiment}_${key}.bam" --write-index
    """
}

process dorado_mod {
    storeDir "bam/sorted_mod"

    label 'gpu'
    label 'std_conda'
    // conda 'bioconda::samtools==1.21'

    input:
        tuple val(experiment), val(key), val(dorado), path(base_model), path(mod_model)

    output:
        path "${experiment}_${key}.bam"
        path "${experiment}_${key}.bam.csi"

    script:
    reference = fetch_ref(experiment)
    minQscore = params.toolConfig.dorado.minQscore
    """
        ${dorado}/bin/dorado basecaller \
            ${base_model} ${params.pod5dir}/${experiment}_5kHz \
            ${params.toolConfig.dorado.run_flags} --reference ${reference} \
            --modified-bases-models ${mod_model} | samtools sort -O BAM -@ ${task.cpus} -o "${experiment}_${key}.bam" --write-index
    """
}
