include { 
    fetch_ref;
    fetch_pod5_loc
} from '../lib/utils.groovy'


process bam_fn_reorder {
    label 'cpu'
    storeDir "bam/sorted_move_fnord"
    when: params.runControls.other_tools.deepbam.size()>0 || params.runControls.other_tools.deeppnat.size()>0

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

    input:  path bam
    output: 
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepbam.tsv"

    script:
    pod5=fetch_pod5_loc("$bam")
    reference=fetch_ref("$bam")
    """
        ${params.toolConfig.deepbam.install_dir} extract_and_call_mods \
            $pod5 $bam $reference DNA ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepbam.tsv \
            ${params.toolConfig.deepbam.model} \
            ${params.toolConfig.deepbam.call_flags}
    """
}

process deepplant_call {
    storeDir "tool_out/deepplant/readwise"
    label 'gpu'

    input:  path bam
    output:
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_cpg.tsv", emit: 'cpg'
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chg.tsv", emit: 'chg'
        path "${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chh.tsv", emit: 'chh'

    script:
    pod5=fetch_pod5_loc("$bam")
    reference=fetch_ref("$bam")
    """
        ${params.toolConfig.deeplant.install_dir} extract_and_call_mods \
            $pod5 $bam $reference DNA ./ \
            ${params.toolConfig.deeplant.model} \
            ${params.toolConfig.deeplant.call_flags}

        mv cpg_result.txt ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_cpg.tsv
        mv chg_result.txt ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chg.tsv
        mv chh_result.txt ${bam.baseName.split('\\.fnTagOrdered')[0]}.deepplant_chh.tsv
    """
}

def which_tool(file, idx){
    return file.split("\\.")[idx]
}

process deetool_aggregate {
    label 'std_conda'

    publishDir { "tool_out/${which_tool("$readwise_file.baseName", -1)}/aggregated" }
    
    label 'cpu'
    label 'std_conda'
    
    input:  path readwise_file
    output: path "${readwise_file.baseName}.aggregated.tsv"
    script:
    """
        python ${projectDir}/scripts/deepbam_aggregate.py --file_path  $readwise_file \
            --aggregation_output ${readwise_file.baseName}.aggregated.tsv \
            --threshold 0.5
    """
}

process deeptool_rebed{
    publishDir { "tool_out/${which_tool("$input_file.baseName", -2)}/rebed" }

    label 'cpu'
    conda 'bioconda::bedtools==2.30.0'

    input:  path input_file
    output: path "${input_file.baseName}.rebed.tsv"

    script:
    reference=fetch_ref("${input_file}")
    """
        [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus} ${reference}
        [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome

        tmp=\$(mktemp /tmp/deeptools_bed_fasta.XXXX);
        bedtools getfasta -fi $reference -bed $input_file -bedOut | \
            awk 'BEGIN{OFS="\\t"} { if(toupper(\$7)=="C") \$7="+"; \
            else if(toupper(\$7)=="G") \$7="-"; print \$1,\$2,\$3,\$5+\$6,".",\$7,\$5,\$6,\$4 }' > \$tmp;
        
        bedtools slop -g ${reference}.genome -l 5 -r 11 -i \$tmp -s | bedtools getfasta -fi ${reference} -bed - -tab -s | cut -f 2 | paste \$tmp - > ${input_file.baseName}.rebed.tsv
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
        python ${projectDir}/scripts/deeptools_consolidate.py $input_file ${input_file.baseName}.std.bed
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
        python ${projectDir}/scripts/deeptools_consolidate.py $input_file ${input_file.baseName}.std.bed
    """
}
