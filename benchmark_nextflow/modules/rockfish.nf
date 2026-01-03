include {
    fetch_ref ;
    fetch_pod5_loc
} from '../lib/utils.groovy'

process setup_rockfish {
    storeDir "${params.tooling_dir}"

    output:
    tuple path("${params.toolConfig.rockfish.install_dir}"), path("${params.toolConfig.rockfish.model_dir}")

    shell:
    """
        source "\$(conda info --base)/etc/profile.d/conda.sh"

        conda create --name "${params.toolConfig.rockfish.conda}" python=3.9 -f -y
        conda activate "${params.toolConfig.rockfish.conda}"
        git clone -b r10.4.1 https://github.com/lbcb-sci/rockfish.git --single-branch "${params.toolConfig.rockfish.install_dir}"
        pip install --extra-index-url https://download.pytorch.org/whl/cu118 ./${params.toolConfig.rockfish.install_dir}/
        
        [ ! -d "${params.toolConfig.rockfish.model_dir}" ] && mkdir -p "${params.toolConfig.rockfish.model_dir}"
        rockfish download -m 5kHz -s "${params.toolConfig.rockfish.model_dir}"
    """
}

process rockfish_call {
    storeDir "rockfish"
    label 'gpu'
    label 'rockfish'

    conda "${params.conda_dir}/envs/${params.toolConfig.rockfish.conda}"

    input:
    tuple path(rockfish), path(rockfishmodels), path(input_bam), path(input_bam_idx)

    output:
    path "${input_bam.baseName}.tsv"

    shell:
    pod5 = fetch_pod5_loc("${input_bam}")
    """
        rockfish inference \
            -i ${pod5} --bam_path ${input_bam} \
            --model_path "${rockfishmodels}/${params.toolConfig.rockfish.model}" \
            -t ${task.cpus} ${params.toolConfig.rockfish.call_flags};
        mv predictions.tsv ${input_bam.baseName}.tsv
    """
}

process rockfish_map_generate {
    storeDir "rockfish/readwise"
    
    label 'cpu'
    label 'rockfish'

    input:
    tuple path(input_bam), path(input_bam_idx), path(prediction_file)

    output:
    path "${input_bam.baseName}.bam_ref_map.tsv"
    path prediction_file

    shell:
    reference = fetch_ref("${input_bam}")
    """
        python ${projectDir}/scripts/rockfish_extract_ref_pos.py --workers ${task.cpus} ${input_bam} ${reference} > ${input_bam.baseName}.bam_ref_map.tsv
    """
}

process rockfish_intersect {
    publishDir "rockfish", mode: "copy"

    label 'cpu'
    label 'std_conda'

    input:
    path ref_map
    path prediction_file

    output:
    path "${prediction_file.baseName}.intercept.tsv"

    script:
    """
        python ${projectDir}/scripts/rockfish_intersect.py -i ${prediction_file} -r ${ref_map} -o "${prediction_file.baseName}.intercept.tsv"
    """
}

process rockfish_aggregate {
    publishDir "rockfish", mode: "copy"

    label 'cpu'
    label 'std_conda'

    input:
    path read_file_mapped

    output:
    path "${read_file_mapped.baseName}.aggregated.tsv"

    script:
    """
        python ${projectDir}/scripts/rockfish_aggregate.py -i ${read_file_mapped} -o "${read_file_mapped.baseName}.aggregated.tsv"
    """
}

process rockfish_getfasta {
    publishDir "rockfish", mode: "copy"

    label 'cpu'
    conda 'bioconda::bedtools==2.30.0'

    input:
    path infile

    output:
    path "${infile.baseName}.rebed.tsv"

    script:
    reference = fetch_ref("${infile}")
    """
        [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus} ${reference}
        [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome

        tmp=\$(mktemp /tmp/rockfish_bed_fasta.XXXX);
        bedtools getfasta -fi ${reference} -bed ${infile} -bedOut | \
            awk 'BEGIN{OFS="\\t"} { if(toupper(\$7)=="C") \$7="+"; else if(toupper(\$7)=="G") \$7="-"; print \$1,\$2,\$3,\$5+\$6,".",\$7,\$5,\$6,\$4 }' > \$tmp;
        bedtools slop -g ${reference}.genome -l 5 -r 11 -s -i \$tmp | \
            bedtools getfasta -s -fi ${reference} -bed - -tab | cut -f 2 | paste \$tmp - > "${infile.baseName}.rebed.tsv";
        
        rm \$tmp;
    """
}

process consolidate_rockfish {
    publishDir "meta/rockfish", mode: "copy"
    label 'cpu'
    label 'std_conda'

    input:
    path input_file

    output:
    path "${input_file.baseName}.std.bed"

    script:
    """
        python ${projectDir}/scripts/rockfish_consolidate.py ${input_file} ${input_file.baseName}.std.bed
    """
}
