include { fetch_ref ; modkitParams } from '../lib/utils.groovy'

process install_modkit {
    storeDir "${params.tooling_dir}"

    output: path "${params.toolConfig.modkit.install_dir}"

    script:
    version=params.toolConfig.modkit.version
    versionurl=params.toolConfig.modkit.version.replace(/-/, '')
    // containerised =  workflow.containerEngine!=null ? "containerised" : ""
    """
        wget https://github.com/nanoporetech/modkit/releases/download/v${version}/modkit_v${versionurl}_u16_x86_64.tar.gz
            tar -xvzf modkit_v${versionurl}_u16_x86_64.tar.gz
            mv dist* ${params.toolConfig.modkit.install_dir}
    """
}

// def fetch_exec(tool, def_exec){
//     def install_dir = params.toolConfig[tool].install_dir
//     install_dir = (install_dir=="" ||  install_dir==null) ? "" : "${params.toolConfig[tool].install_dir}/"

//     def exec = params.toolConfig[tool].executable
//     exec = (exec=="" || exec==null) ? def_exec : params.toolConfig[tool].executable

//     return "${install_dir}${exec}"
// }

process modkit_pileup {
    publishDir 'tool_out/dorado/pileup'
    label 'cpu'

    input:
        tuple val(modkit),
            path(input_file), 
            path(input_index)

    output:
        path "${input_file.baseName}.bed"

    script:
    reference = fetch_ref("${input_file}")
    modParams = modkitParams(input_file)
    flags = params.toolConfig.modkit.call_flags==null ? "" : params.toolConfig.modkit.call_flags
    exec  = workflow.containerEngine!=null ? 'modkit' : file("${modkit}/${ params.toolConfig.modkit.executable}")
    """
        ${exec} pileup ${input_file} ${input_file.baseName}.bed -t ${task.cpus} --ref ${reference} ${modParams} ${flags}
    """
}

process modkit_add_ref {
    publishDir 'tool_out/dorado/ref'

    label 'cpu'

    input:
        tuple path(input_file), path(reference), path(faidx), path(genome)

    output:
        path "${input_file.baseName}.ref.bed"

    script:
    """
        bedtools slop -s -l 5 -r 11 -g ${genome} -i ${input_file} \
        | bedtools getfasta -fi ${reference} -bed - -tab -s \
        | cut -f 2 | paste -d '\\t' ${input_file} -  > ${input_file.baseName}.ref.bed
    """
}

process standardise_dorado {
    publishDir "meta/dorado"

    label 'cpu'
    label 'std_conda'

    input:
        path input_file

    output:
        path "${input_file.baseName}.std.bed"

    script:
    """
        modkit_consolidate.py ${input_file} ${input_file.baseName}.std.bed
    """
}
