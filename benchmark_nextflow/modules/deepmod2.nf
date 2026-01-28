include {
    fetch_ref ;
    fetch_deepmod2_model ;
    fetch_pod5_loc
} from '../lib/utils.groovy'

process setup_deepmod2 {
    storeDir "${params.tooling_dir}"

    label 'install'
    
    output:
    path "${params.toolConfig.deepmod2.install_dir}"

    shell:
    """
        source "\$(conda info --base)/etc/profile.d/conda.sh"

        git clone https://github.com/WGLab/DeepMod2.git "${params.toolConfig.deepmod2.install_dir}"

        conda env create \
            -f "${params.toolConfig.deepmod2.install_dir}"/environment.yml \
            -n "${params.toolConfig.deepmod2.conda}" -y

        ## installing cuda packages

        conda activate "${params.toolConfig.deepmod2.conda}"
        cuda_version=\$(nvidia-smi | grep -oP 'CUDA Version: \\K[\\d.]+' | sed 's/\\.//')
        echo installing \$cuda_version
        pip install torch torchvision --index-url https://download.pytorch.org/whl/cu\$cuda_version
    """
}

process deepmod2_call {
    storeDir "tool_out/deepmod2"

    label 'gpu' 
    label 'deepmod2'

    input:
    tuple path(deepmod2_dir), path(input_bam), path(input_bam_idx), path(reference), val(model)

    output:
    path "${input_bam.baseName}", emit: folder

    script:
    model_path = fetch_deepmod2_model("${input_bam}", "${model}")
    reference = fetch_ref("${input_bam}")
    pod5 = fetch_pod5_loc("${input_bam}")
    """
        python ${deepmod2_dir}/deepmod2 detect \
            --bam ${input_bam} \
            --input ${pod5} \
            --model ${model_path} \
            --file_type pod5 \
            --threads ${task.cpus} \
            --ref ${reference} \
            --seq_type dna \
            --output ${input_bam.baseName} \
            ${params.toolConfig.deepmod2.call_flags}
    """
}

process deepmod2_add_ref {
    publishDir "tool_out/deepmod2", mode: "copy"

    label 'cpu'

    input:
    tuple path(input_folder), val(model)

    output:
    tuple path("${input_folder}_${model}.deepmod2.aggregated.rebed.ref.tsv"), val(model)

    script:
    reference = fetch_ref("${input_folder}")
    """
        [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus} ${reference}
        [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome

        tmp=\$(mktemp /tmp/deepmod2_metadata.XXXX);
        sed '1d' ${input_folder}/output.per_site | awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$6,\$5,\$4,\$7,\$8,\$9}' > \$tmp;

        bedtools slop -g ${reference}.genome -l 5 -r 11 -i \$tmp -s | bedtools getfasta -fi ${reference} -bed - -tab -s | cut -f 2 | paste \$tmp - > "${input_folder}_${model}.deepmod2.aggregated.rebed.ref.tsv"
        rm \$tmp
    """
}

process consolidate_deepmod2 {
    publishDir "meta/deepmod2", mode: "copy"

    label 'cpu'
    label 'std_conda'

    input:
    tuple path(input_file), val(model)

    output:
    path "${input_file.baseName}.std.bed"

    script:
    """
        python ${projectDir}/scripts_common/deepmod2_consolidate.py ${input_file} ${input_file.baseName}.std.bed ${model}
    """
}
