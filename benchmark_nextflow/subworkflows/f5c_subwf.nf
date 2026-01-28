include { pod5_to_blow5 ; f5c_idx_fastq ; f5c_call ; f5c_restrand } from '../modules/f5c.nf'
include { setup_f5c ; f5c_aggregate ; f5c_add_fasta ; consolidate_f5c } from '../modules/f5c.nf'
include {
    f5c_aggregate as f5c_aggregate_stranded ;
    f5c_add_fasta as f5c_add_fasta_stranded ;
    consolidate_f5c as consolidate_f5c_stranded;
} from '../modules/f5c.nf'

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

        setup_f5c |                           // install f5c if not exist
                combine(pod5_to_blow5.output) |   // map on blow5 folders to fastq
                f5c_idx_fastq                     // index fastq and blow5

        f5c_idx_fastq.output                  // map blow5 folders and call
            .map{ it -> ["${it[3].toString().split('/')[-1].split('_5kHz')[0]}_5kHz", *it ]}
            .combine(reference_map_ob, by: 0)
            .map{ it -> [ *it[1..9], file(it[-1]) ] } | f5c_call

    emit:
        setup_f5c = setup_f5c.output
        f5c_call = f5c_call.output
}

workflow F5C_UNSTRANDED_WORKFLOW {
    take:
        setup_f5c_ch
        f5c_call_ch

    main:
        f5c_call_ob  = setup_f5c_ch.combine(f5c_call_ch)
        f5c_aggregate(f5c_call_ob)
        f5c_add_fasta(f5c_aggregate.output)
        consolidate_f5c(f5c_add_fasta.output)
    // emit:
    //     f5c_restrand = f5c_restrand.output
}

workflow F5C_RESTRANDED_WORKFLOW {
    take:
        setup_f5c_ch
        f5c_call_ch

    main:
        f5c_restrand(f5c_call_ch)
        f5c_call_ob  = setup_f5c_ch.combine(f5c_restrand.output)
        f5c_aggregate_stranded(f5c_call_ob)
        f5c_add_fasta_stranded(f5c_aggregate_stranded.output)
        consolidate_f5c_stranded(f5c_add_fasta_stranded.output)

    // emit:
    //     f5c_restrand = f5c_restrand.output
}
