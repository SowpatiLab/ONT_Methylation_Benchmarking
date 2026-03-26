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

    script:
    """
        git clone https://github.com/nanoporetech/rerio.git rerio
    """
}

process download_rerio_models {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.model_dir}"

    input:  tuple val(base_model), val(key),  val(model_name), val(dorado), val(rerio)
    output: tuple val(base_model), val(key), path(model_name)

    script:
    containerised = workflow.containerEngine!=null ? "containerised" : "conda"
    if(workflow.containerEngine!=null) {println("containerised in ${workflow.containerEngine} | ${containerised}")}
    """
        if [[ $containerised == "containerised" ]];
        then
            echo $containerised
            echo 'containerised rerio download';
            cp -r /tooling/rerio/download_model.py /tooling/rerio/dorado_models . ;
            ./download_model.py --dorado dorado_models/${model_name}_url ;
            mv dorado_models/${model_name} . ;
            rm -r download_model.py dorado_models ;
        else
            echo 'conda rerio download';
            ${rerio}/download_model.py --dorado ${rerio}/dorado_models/${model_name}_url
            cp -r ${rerio}/dorado_models/${model_name} .
        fi
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
    dorver="${params.toolConfig.dorado.dorado_version_aliases[dorado_version_map]}"
    """
        [ ! -f ${dorado_version_map}-linux-x64.tar.gz ] && wget https://cdn.oxfordnanoportal.com/software/analysis/${dorado_version_map}-linux-x64.tar.gz \
            -O ${dorado_version_map}-linux-x64.tar.gz
        tar -xvzf ${dorado_version_map}-linux-x64.tar.gz 
        mv "${dorado_version_map}"-linux-x64 "${dorver}"
    """
}

process download_dorado_base_model {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.model_dir}"

    input:
        tuple val(model_name), val(dorado)

    output:
        tuple val(model_name), val(dorado), path(model_name)

    script:
    dorado_executabel = workflow.containerEngine!=null ? "${dorado}" : "${dorado}/bin/dorado"
    """
        ${dorado_executabel} download --model ${model_name};
    """
}

process download_dorado_mod_models {
    label 'download'
    storeDir "${params.tooling_dir}/${params.toolConfig.dorado.model_dir}"

    input:
        tuple val(base_model), val(key), val(model_name), val(dorado)

    output:
        tuple val(base_model), val(key), path(model_name)

    script:
    dorado_executabel = workflow.containerEngine!=null ? "${dorado}" : "${dorado}/bin/dorado"
    """
        ${dorado_executabel} download --model ${model_name};
    """
}

process dorado_move {
    publishDir "${params.output_dir}/bam/sorted_mod"
    
    label 'std_conda'
    label 'gpu'

    input:
        tuple val(experiment), val(key), val(dorado), path(base_model), path(pod5), path(reference)

    output:
        tuple path("${experiment}_${key}.bam"), path("${experiment}_${key}.bam.csi")

    script:
    dorado_executabel = workflow.containerEngine!=null ? "${dorado}" : "${dorado}/bin/dorado"
    """
        ${dorado_executabel} basecaller \
            ${base_model} ${pod5} \
            ${params.toolConfig.dorado.run_flags} \
            --reference ${reference} \
            --emit-moves \
            | samtools sort -O BAM -@ ${task.cpus} -o "${experiment}_${key}.bam" --write-index;
    """
}

process dorado_mod {
    publishDir "${params.output_dir}/bam/sorted_mod"

    label 'std_conda'
    label 'gpu'
    
    input:
        tuple val(experiment), val(key), val(dorado), path(base_model), path(mod_model), path(pod5), path(reference)

    output:
        tuple path("${experiment}_${key}.bam"), path("${experiment}_${key}.bam.csi")

    script:
    dorado_executabel = workflow.containerEngine!=null ? "${dorado}" : "${dorado}/bin/dorado"
    """
        ${dorado_executabel} basecaller \
            ${base_model} ${pod5} \
            ${params.toolConfig.dorado.run_flags} \
            --modified-bases-models ${mod_model} \
            --reference ${reference} \
            | samtools sort -O BAM -@ ${task.cpus} -o "${experiment}_${key}.bam" --write-index;
    """
}