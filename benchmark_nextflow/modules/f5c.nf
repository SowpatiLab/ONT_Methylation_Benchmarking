include {
    fetch_ref ;
    get_sample_name ;
    fetch_pod5_loc
} from '../lib/utils.groovy'

process bam_to_fastq {
    publishDir "intermediary/fastq", mode: "copy", pattern: "*.fastq*"

    label 'cpu'
    label 'std_conda'

    input:
        tuple path(bam), path(bam_idx)

    output:
        tuple path(bam), path(bam_idx), path("${bam.baseName}.fastq")

    script:
    """
        samtools fastq -@ ${task.cpus} ${bam} -0 ${bam.baseName}.fastq
    """
}

process pod5_to_blow5 {
    publishDir "intermediary/blow5", mode: "copy", pattern: "*.blow5*"

    label 'cpu'
    label 'std_conda'

    input:
        tuple path(pod5), path(bam), path(bam_idx), path(fastq), val(key)

    output:
        tuple path("*.blow5"), path(bam), path(bam_idx), path(fastq)

    script:
    """
        blue-crab p2s ${pod5} -t ${task.cpus} -o ${key}.blow5
    """
}

process setup_f5c {
    storeDir "${params.tooling_dir}"

    output:
    path "${params.toolConfig.f5c.install_dir}"

    shell:
    """
        [ ! -d ${params.toolConfig.f5c.install_dir} ] && mkdir ${params.toolConfig.f5c.install_dir}
        wget "https://github.com/hasindu2008/f5c/releases/download/${params.toolConfig.f5c.version}/f5c-${params.toolConfig.f5c.version}-binaries.tar.gz" && \
            tar xvf f5c-${params.toolConfig.f5c.version}-binaries.tar.gz -C ${params.toolConfig.f5c.install_dir} --strip-component=1
    """
}

process f5c_idx_fastq {
    publishDir "intermediary/fastq", mode: "copy", pattern: "*.fastq*"
    publishDir "intermediary/blow5", mode: "copy", pattern: "*.blow5*"

    input:
    tuple path(install_dir), path(blow5), path(bam), path(bam_idx), path(fastq)

    output:
    tuple path(install_dir), path(blow5), path("${blow5}.idx"), path(bam), path(bam_idx), path(fastq), path("${fastq}.index"), path("${fastq}.index.fai"), path("${fastq}.index.gzi")

    script:
    """
        ${install_dir}/f5c_x86_64_linux_cuda index --slow5 ${blow5} ${fastq};
    """
}

process f5c_call {
    storeDir "tool_out/f5c/readwise"
    label 'gpu'
    
    input:
    tuple path(install_dir), path(blow5), path(blow5_idx), path(bam), path(bam_idx), path(fastq), path(fastq_index), path(fastq_index_fai), path(fastq_index_gzi), path(reference)

    output:
    path "${bam.baseName}.f5c.tsv"

    script:
    execut = "${params.toolConfig.f5c.executable}"
    """
        ${install_dir}/${execut} call-methylation \
            ${params.toolConfig.f5c.call_flags} \
            -t ${task.cpus} \
            --slow5 ${blow5} \
            -b ${bam} \
            -g ${reference} \
            -r ${fastq} > ${bam.baseName}.f5c.tsv
    """
}

process f5c_restrand {
    publishDir "tool_out/f5c/readwise_stranded", mode: "copy", pattern: "*.tsv"
    
    label 'cpu'
    label 'std_conda'

    input:
    path input_file

    output:
    path "${input_file.baseName}.stranded.tsv"

    when:
    params.runControls.other_tools.f5c_stranded.size() > 0

    script:
    reference = fetch_ref("${input_file}")
    """
        python ${projectDir}/scripts/f5c_restrand.py -i ${input_file} -r ${reference} -o ${input_file.baseName}.stranded.tsv
    """
}

// todo: batch wise merging
// process merge_f5c_calls {
//     input:  tuple path(install_dir), path(input_files)
//     output: tuple path(install_dir), path("merged.f5c.tsv")
//     script:
//     """
//         python ${projectDir}/scripts/merge_f5c_calls.py -i ${input_files.join(' ')} -o merged.f5c.tsv
//     """
// }

def isStranded(filename) {
    return filename.toString().contains("stranded")
}

process f5c_aggregate {
    publishDir { "tool_out/f5c/aggregated${isStranded(input_file) ? '_stranded' : ''}" }
    
    label 'cpu'
    label 'std_conda'

    input:
    tuple path(install_dir), path(input_file)

    output:
    path "${input_file.baseName}.aggregated.tsv"

    script:
    reference = fetch_ref("${input_file}")
    execut = "${params.toolConfig.f5c.executable}"
    """
        ${install_dir}/${execut}  meth-freq -i ${input_file} -s > ${input_file.baseName}.aggregated.tsv
    """
}

process f5c_add_fasta {
    publishDir { "tool_out/f5c/aggregated_fasta${isStranded(input_file) ? '_stranded' : ''}" }

    label 'cpu'
    conda 'bioconda::bedtools==2.30.0'

    input:
    path input_file

    output:
    path "${input_file.baseName}.rebed.ref.tsv"

    script:
    reference = fetch_ref("${input_file}")
    stranded = isStranded("${input_file}")
    if (stranded) {
        """
            [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus} ${reference}
            [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome

            echo "running started pipeline"
            tmp=\$(mktemp /tmp/f5c_bed_fasta_stranded.XXXX);
                    
            bedtools getfasta -fi ${reference} -bed <(tail -n +2 ${input_file}) -bedOut | \
                awk 'BEGIN{OFS="\\t"} { if(toupper(\$9)=="C") \$9="+"; else if(toupper(\$9)=="G") \$9="-"; print \$1,\$2,\$3,\$5,".",\$9,\$6,\$5-\$6,\$7 }' > \$tmp;

            bedtools slop -l 5 -r 11 -s -g ${reference}.genome -i \$tmp \
                | bedtools getfasta -fi ${reference} -bed - -tab -s | cut -f 2 | paste \$tmp - > ${input_file.baseName}.rebed.ref.tsv;
            rm \$tmp
        """
    }
    else {
        """
            [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus} ${reference}
            [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome

            echo "running un-started pipeline"
            tmp=\$(mktemp /tmp/f5c_bed_fasta.XXXX);

            sed '1d' ${input_file} | awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$2+1,\$5,".","+",\$6,\$5-\$6,\$7}' > \$tmp;
            bedtools slop -l 5 -r 11 -s -g ${reference}.genome -i \$tmp \
                | bedtools getfasta -fi ${reference} -bed - -tab -s | cut -f 2 | paste \$tmp - > ${input_file.baseName}.rebed.ref.tsv;
            rm \$tmp
        """
    }
}

process consolidate_f5c {
    publishDir { "meta/f5c${isStranded(input_file) ? '_stranded' : ''}" }, mode: "copy"

    label 'cpu'
    label 'std_conda'

    input:
    path input_file

    output:
    path "${input_file.baseName}.std.bed"

    script:
    """
        python ${projectDir}/scripts/f5c_consolidate.py ${input_file} ${input_file.baseName}.rebed.ref.bed
    """
}
