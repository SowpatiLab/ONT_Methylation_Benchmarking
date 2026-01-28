include { 
    fetch_ref;
} from '../lib/utils.groovy'


process bam_fn_reorder {
    label 'cpu'
    storeDir "bam/sorted_move_fnord"

    input:  tuple path(bam), path(bam_idx)
    output: path "${bam.baseName}.fnTagOrdered.bam"
    
    shell:
    """
        samtools sort -t fn $bam -@ ${task.cpus} -O BAM -o "${bam.baseName}.fnTagOrdered.bam"
    """
}

process deepbam_call {
    storeDir "tool_out/deepbam/readwise"
    label 'gpu'
    conda "${params.conda_dir}/envs/${params.toolConfig.deepbam.conda}"

    container '/data1/ccmb/reuben/benchmarking/bacterial_nf/apptainer_final/ontMethylationBenchmarking.sif'
    containerOptions '--nv'

    input: tuple path(pod5), path(bam), path(reference), path(model)
    output: 
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepbam.tsv"

    script:
    executable = workflow.containerEngine!=null ? "DeepBAM" : "${params.toolConfig.deepbam.install_dir}/${params.toolConfig.deepbam.executable}"
    model_name = workflow.containerEngine!=null ? "/DeepBAM/traced_script_module/${model}" :  model
    println("running deepbam_call on: ${workflow.containerEngine=="apptainer" ? workflow.containerEngine : task.executor}")

    """
        ${executable} extract_and_call_mods \
            ${pod5} ${bam} ${reference} \
            DNA ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepbam.tsv \
            ${model_name} \
            ${params.toolConfig.deepbam.call_flags}
    """
}

process deepplant_call {
    storeDir "tool_out/deepplant/readwise"
    label 'gpu'

    container '/data1/ccmb/reuben/benchmarking/bacterial_nf/apptainer_final/ontMethylationBenchmarking.sif'
    containerOptions '--nv'

    input:  tuple path(pod5), path(bam), path(reference), path(model)
    output:
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_cpg.tsv", emit: 'cpg'
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chg.tsv", emit: 'chg'
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chh.tsv", emit: 'chh'

    script:
    executable = workflow.containerEngine!=null ? "DeepPlant" : "${params.toolConfig.deepplant.install_dir}/${params.toolConfig.deepplant.executable}"
    model_name = workflow.containerEngine!=null ? "/DeepPlant/model/bilstm/${model}" :  model
    println("running deepplant_call on: ${workflow.containerEngine=="apptainer" ? workflow.containerEngine : task.executor}")
    """
        ${executable} extract_and_call_mods \
            ${pod5} ${bam} ${reference} DNA ./ ${model} \
            ${params.toolConfig.deepplant.call_flags}

        mv cpg_result.txt ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_cpg.tsv
        mv chg_result.txt ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chg.tsv
        mv chh_result.txt ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chh.tsv
    """
}

def which_tool(file, idx){
    return file.split("\\.")[idx].split('_')[0]
}

process deetool_aggregate {
    label 'std_conda'

    publishDir { "tool_out/${which_tool("${readwise_file.baseName}", -1)}/aggregated" }
    
    label 'cpu'
    label 'std_conda'
    
    input:  path readwise_file
    output: path "${readwise_file.baseName}.aggregated.tsv"
    script:
    aggregation_flags=params.toolConfig[which_tool("${readwise_file.baseName}", -1)].aggregation_flags
    """
        python ${projectDir}/scripts_common/deepbam_aggregate.py --file_path  $readwise_file \
            --aggregation_output ${readwise_file.baseName}.aggregated.tsv \
            $aggregation_flags
    """
}

process deeptool_rebed{
    publishDir { "tool_out/${which_tool("$input_file.baseName", -2)}/rebed" }

    label 'cpu'
    label 'std_conda'
    // conda 'bioconda::bedtools==2.30.0'

    input:  path input_file
    output: path "${input_file.baseName}.rebed.ref.tsv"

    script:
    reference=fetch_ref("${input_file}")
    """
        [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus} ${reference}
        [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome

        tmp=\$(mktemp /tmp/deeptools_bed_fasta.XXXX);
        bedtools getfasta -fi $reference -bed $input_file -bedOut | \
            awk 'BEGIN{OFS="\\t"} { if(toupper(\$7)=="C") \$7="+"; \
            else if(toupper(\$7)=="G") \$7="-"; print \$1,\$2,\$3,\$5+\$6,".",\$7,\$5,\$6,\$4 }' > \$tmp;
        
        bedtools slop -g ${reference}.genome -l 5 -r 11 -i \$tmp -s | bedtools getfasta -fi ${reference} -bed - -tab -s | cut -f 2 | paste \$tmp - > ${input_file.baseName}.rebed.ref.tsv
    """
}

process deepbam_consolidate {
    publishDir "meta/deepbam", mode: "copy"

    label 'cpu'
    label 'std_conda'

    input:  path input_file
    output: path "${input_file.baseName}.std.bed"
    script:
    """
        python ${projectDir}/scripts_common/deeptools_consolidate.py $input_file ${input_file.baseName}.std.bed
    """
}

process deepplant_consolidate {
    publishDir "meta/deepplant", mode: "copy"
    label 'cpu'
    label 'std_conda'
    
    input:  path input_file
    output: path "${input_file.baseName}.std.bed"
    script:
    """
        python ${projectDir}/scripts_common/deeptools_consolidate.py $input_file ${input_file.baseName}.std.bed
    """
}
