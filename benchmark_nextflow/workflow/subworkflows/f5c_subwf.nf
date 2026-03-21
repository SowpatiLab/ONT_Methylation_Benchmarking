include { pod5_to_blow5 ; f5c_idx_fastq ; f5c_call ; f5c_restrand } from '../modules/f5c.nf'
include { setup_f5c ; f5c_aggregate ; f5c_add_fasta ; consolidate_f5c } from '../modules/f5c.nf'
include {
    f5c_aggregate as f5c_aggregate_stranded ;
    f5c_add_fasta as f5c_add_fasta_stranded ;
    consolidate_f5c as consolidate_f5c_stranded;
} from '../modules/f5c.nf'

include { attachReferenceFiles } from '../lib/utils.groovy'

workflow F5C_GENERAL_WORKFLOW {
    take:
        fastq_ch
        reference_map_ob

    main: 
        fastq_ch
            .map{ bam, csi, fq -> ["${bam.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", bam, csi, fq ]}
            .combine(reference_map_ob, by: 0)
            .map{ key, bam, csi, fq, p5, fa -> [p5, bam, csi, fq,  key] } |
            pod5_to_blow5                     // blow5 generation

        f5c_exe = Channel.empty()
        if(workflow.containerEngine!=null){
            f5c_exe = Channel.of("${params.toolConfig.f5c.executable}")
        } else {
            setup_f5c()
            f5c_exe = setup_f5c.output
        }
        f5c_exe
            .combine(pod5_to_blow5.output) |   // map on blow5 folders to fastq
            f5c_idx_fastq                     // index fastq and blow5

        f5c_idx_fastq.output                  // map blow5 folders and call
            .map{ it -> ["${it[3].toString().split('/')[-1].split('_5kHz')[0]}_5kHz", *it ]}
            .combine(reference_map_ob, by: 0)
            .map{ it -> [ *it[1..9], file(it[-1]) ] } | f5c_call

    emit:
        setup_f5c = f5c_exe
        f5c_call = f5c_call.output
}

workflow F5C_UNSTRANDED_WORKFLOW {
    take:
        reference_map_ch
        setup_f5c_ch
        f5c_call_ch

    main:
        f5c_call_ob  = setup_f5c_ch.combine(f5c_call_ch)
        f5c_aggregate(f5c_call_ob)

        attachReferenceFiles(f5c_aggregate.output, reference_map_ch) | 
                f5c_add_fasta |
                consolidate_f5c
    // emit:
    //     f5c_restrand = f5c_restrand.output
}

workflow F5C_RESTRANDED_WORKFLOW {
    take:
        reference_map_ch
        setup_f5c_ch
        f5c_call_ch

    main:
        def reference_map = new groovy.yaml.YamlSlurper().parseText(file(params.references).text)
        f5c_call_ch
            .map{ bed -> [ 
                    bed, file(reference_map_ch['references'][bed.baseName.split('_')[0]]),
            ]} | 
            f5c_restrand
            
        setup_f5c_ch
            .combine(f5c_restrand.output) |
            f5c_aggregate_stranded
        
        attachReferenceFiles(f5c_aggregate_stranded.output, reference_map_ch) | 
            f5c_add_fasta_stranded |
            consolidate_f5c_stranded

    emit:
        f5c_restrand = f5c_restrand.output
}
