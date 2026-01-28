include { fetch_ref ; modkitParams } from '../lib/utils.groovy'

process install_modkit {
    storeDir "${params.tooling_dir}"

    output: path "${params.toolConfig.modkit.install_dir}"
    script:
    version=params.toolConfig.modkit.version
    versionurl=params.toolConfig.modkit.version.replace(/-/, '')
    """
        if ${params.tooling_dir}/modkit &> /dev/null
        then
            exit 0
        else 
            wget https://github.com/nanoporetech/modkit/releases/download/v${version}/modkit_v${versionurl}_u16_x86_64.tar.gz
            tar -xvzf modkit_v${versionurl}_u16_x86_64.tar.gz
            mv dist* ${params.toolConfig.modkit.install_dir}
        fi
    """
}

def fetch_exec(tool, def_exec){
    def install_dir = params.toolConfig[tool].install_dir
    install_dir = (install_dir=="" ||  install_dir==null) ? "" : "${params.toolConfig[tool].install_dir}/"

    def exec = params.toolConfig[tool].executable
    exec = (exec=="" || exec==null) ? def_exec : params.toolConfig[tool].executable

    return "${install_dir}${exec}"
}

process modkit_pileup {
    publishDir 'tool_out/dorado/pileup'
    label 'cpu'

    input:
        tuple path(modkit), path(input_file), path(input_index)

    output:
        path "${input_file.baseName}.bed"

    script:
    reference = fetch_ref("${input_file}")
    modParams = modkitParams(input_file)
    flags = params.toolConfig.modkit.call_flags==null ? "" : params.toolConfig.modkit.call_flags
    exec  = fetch_exec('modkit', 'modkit')
    """
        ${exec} pileup ${input_file} ${input_file.baseName}.bed -t ${task.cpus} --ref ${reference} ${modParams} ${flags}
    """
}

process modkit_add_ref {
    storeDir 'tool_out/dorado/ref'

    label 'cpu'

    input:
        path input_file

    output:
        path "${input_file.baseName}.ref.bed"

    script:
    reference = fetch_ref("${input_file}")
    """
        [ ! -d ${reference}.fai ] && samtools faidx -@ ${task.cpus}
        [ ! -d ${reference}.genome ] && cut -f 1-2 ${reference}.fai > ${reference}.genome
        
        bedtools slop -s -l 5 -r 11 -g ${reference}.genome -i ${input_file} | bedtools getfasta -fi ${reference} -bed - -tab -s | cut -f 2 | paste -d '\t' ${input_file} - > ${input_file.baseName}.ref.bed
    """
}

process standardise_dorado {
    storeDir "meta/dorado"

    label 'cpu'
    label 'std_conda'

    input:
        path input_file

    output:
        path "${input_file.baseName}.std.bed"

    script:
    """
        python ${projectDir}/scripts_common/modkit_consolidate.py ${input_file} ${input_file.baseName}.std.bed
    """
}
