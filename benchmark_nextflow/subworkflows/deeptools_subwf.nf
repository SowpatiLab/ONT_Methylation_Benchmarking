
include { bam_fn_reorder; deepbam_call; deepplant_call } from '../modules/deeptools.nf'
include {
    deetool_aggregate as deepbam_aggregate;
    deeptool_rebed as deepbam_rebed;
    deepbam_consolidate
} from '../modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_cpg; deeptool_rebed as deepplant_rebed_cpg; deepplant_consolidate as deepplant_consolidate_cpg } from '../modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_chg; deeptool_rebed as deepplant_rebed_chg; deepplant_consolidate as deepplant_consolidate_chg } from '../modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_chh; deeptool_rebed as deepplant_rebed_chh; deepplant_consolidate as deepplant_consolidate_chh } from '../modules/deeptools.nf'


workflow DEEPTOOLS_GENERAL_SUBWORKFLOW {
    take:
        samtools_cleansed_ch
    main: 
        bam_fn_reorder(samtools_cleansed_ch)
    emit: 
        bam_fn_reorder = bam_fn_reorder.output
}

workflow DEEPBAM_SUBWORKFLOW {
    take:
        bam_fn_reorder_ch
        reference_map_ob
        model_ch
    main:
        deepBamIn = bam_fn_reorder_ch
            .map{ it -> [
                "${it.toString().split('/')[-1].split('_5kHz')[0]}_5kHz",
                it
            ] }
            .combine(reference_map_ob, by: 0)
            .combine(model_ch)
            .map{ key, bam, pod5, fasta, model -> [ pod5, bam, fasta, file(model) ] }
                    
        deepbam_call(deepBamIn) 
        deepbam_aggregate(deepbam_call.output)
        deepbam_rebed(deepbam_aggregate.output)
        deepbam_consolidate(deepbam_rebed.output)
}

workflow DEEPPLANT_SUBWORKFLOW {
    take:
        bam_fn_reorder_ch
        reference_map_ob
        model_ch
    main:
        deepPlantIn = bam_fn_reorder_ch
            .map{ it -> [ "${it.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", it ] }
            .combine(reference_map_ob, by: 0)
            .combine(model_ch)
            .map{ key, bam, pod5, fasta, model -> [ pod5, bam, fasta, file(model) ] }

        deepplant_call(deepPlantIn)
        
        deepplant_aggregate_cpg(deepplant_call.output.cpg)
        deepplant_aggregate_chg(deepplant_call.output.chg)
        deepplant_aggregate_chh(deepplant_call.output.chh)

        deepplant_rebed_cpg(deepplant_aggregate_cpg.output)
        deepplant_rebed_chg(deepplant_aggregate_chg.output)
        deepplant_rebed_chh(deepplant_aggregate_chh.output)

        deepplant_consolidate_cpg(deepplant_rebed_cpg.output)
        deepplant_consolidate_chg(deepplant_rebed_chg.output)
        deepplant_consolidate_chh(deepplant_rebed_chh.output)
}