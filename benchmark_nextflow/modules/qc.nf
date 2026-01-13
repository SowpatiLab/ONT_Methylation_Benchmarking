process nanostat {
    publishDir "qc/nanostat", mode: "copy"
    conda "bioconda::nanostat==1.6.0"
    label "cpu"
    label 'std_conda'

    input: tuple path(input_bam), path(input_bam_idx)
    output: 
        path "${input_bam.baseName}.nanostat"
    
    shell:
    """
        NanoStat -t ${task.cpus} --bam ${input_bam} > ${input_bam.baseName}.nanostat
    """
}

process nanoplot {
    publishDir "qc/nanoplot", mode: "copy"

    label "cpu"
    label 'std_conda'
    // conda "bioconda::nanoplot==1.46.0"

    errorStrategy 'retry'
    maxRetries 3

    input: tuple path(input_bam), path(input_bam_idx)
    output: 
        path "${input_bam.baseName}"
    
    shell:
    """
        NanoPlot -t ${task.cpus} --no_static --only-report --minqual 0 --bam ${input_bam} --raw --tsv_stats -o ${input_bam.baseName}
    """
}

// process nanoq {
//     publishDir "qc/nanoq", mode: "copy"

//     label "cpu"
//     label 'std_conda'
//     //conda "bioconda::nanoq==0.10.0"

//     input: tuple path(input_fastq)
//     output: 
//         path "${input_fastq.baseName}.nanoq"
    
//     shell:
//     """
//         nanoq -s -i ${input_fastq} > ${input_fastq.baseName}.nanoq
//     """
// }

process mosdepth {
    publishDir "qc/mosdepth", mode: "copy"
    conda "bioconda::mosdepth=0.3.12"
    label "cpu"

    input: tuple path(input_bam), path(input_bam_idx)
    output: 
        path "${input_bam.baseName}"
    
    shell:
    """
        mkdir ${input_bam.baseName}
        mosdepth -t ${task.cpus} ${input_bam.baseName}/${input_bam.baseName} ${input_bam}
    """
}
