include {
    fetchStdChromosomes
} from '../lib/utils.groovy'

process samtools_move_cleanse {
    storeDir "bam/sorted_move_cleansed"

    label 'std_conda'
    // conda 'bioconda::samtools==1.21'

    input:
        tuple path(input_file), path(input_index)

    output:
        tuple path("${input_file.baseName}.cleansed.bam"), path("${input_file.baseName}.cleansed.bam.csi")

    script:
    chroms = fetchStdChromosomes(input_file)
    """
        samtools view ${input_file} ${chroms} -hbo ${input_file.baseName}.cleansed.bam -@ ${task.cpus} ${params.toolConfig.samtools.cleanse_flags} --write-index
    """
}

process samtools_mod_cleanse {
    storeDir "bam/sorted_mod_cleansed"

    label 'std_conda'

    input:
    path input_file
    path input_index

    output:
    tuple path("${input_file.baseName}.cleansed.bam"), path("${input_file.baseName}.cleansed.bam.csi")

    script:
    chroms = fetchStdChromosomes(input_file)
    """
        samtools view ${input_file} ${chroms} -hbo ${input_file.baseName}.cleansed.bam -@ ${task.cpus} ${params.toolConfig.samtools.cleanse_flags} --write-index
    """
}
