
include { bam_fn_reorder; deepbam_call; deepplant_call } from '../modules/deeptools.nf'
include {
    deetool_aggregate as deepbam_aggregate;
    deeptool_rebed as deepbam_rebed;
    deepbam_consolidate;
    download_deepbam_models;
    download_deepplant_models;
} from '../modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_cpg; deeptool_rebed as deepplant_rebed_cpg; deepplant_consolidate as deepplant_consolidate_cpg } from '../modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_chg; deeptool_rebed as deepplant_rebed_chg; deepplant_consolidate as deepplant_consolidate_chg } from '../modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_chh; deeptool_rebed as deepplant_rebed_chh; deepplant_consolidate as deepplant_consolidate_chh } from '../modules/deeptools.nf'

include { attachReferenceFiles } from '../lib/utils.groovy'

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
        reference_map_ch
        bam_fn_reorder_ch
        reference_map_ob
    main:

        download_deepbam_models()
        deepbam_model = download_deepbam_models
            .output
            .map { f -> file("${f}/${params.toolConfig.deepbam.model}")}

        deebam_in_ch = bam_fn_reorder_ch
            .map{ it -> [
                "${it.toString().split('/')[-1].split('_5kHz')[0]}_5kHz",
                it
            ] }
            .combine(reference_map_ob, by: 0)
            .combine(deepbam_model)
            .map{ key, bam, pod5, fasta, model -> [ pod5, bam, fasta, model ] } | 
            deepbam_call | 
            deepbam_aggregate

        attachReferenceFiles(deepbam_aggregate.output, reference_map_ch) |
            deepbam_rebed |
            deepbam_consolidate 
}

workflow DEEPPLANT_SUBWORKFLOW {
    take:
        reference_map_ch
        bam_fn_reorder_ch
        reference_map_ob
    main:

        download_deepplant_models()
        model_ch = download_deepplant_models
            .output
            .map { f -> file("${f}/${params.toolConfig.deepplant.model}")}

        deepPlantIn = bam_fn_reorder_ch
            .map{ it -> [ "${it.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", it ] }
            .combine(reference_map_ob, by: 0)
            .combine(model_ch)
            .map{ key, bam, pod5, fasta, model -> [ pod5, bam, fasta, file(model) ] }

        deepplant_call(deepPlantIn)
        
        deepplant_aggregate_cpg(deepplant_call.output.cpg)
        deepplant_aggregate_chg(deepplant_call.output.chg)
        deepplant_aggregate_chh(deepplant_call.output.chh)

        attachReferenceFiles(deepplant_aggregate_cpg.output, reference_map_ch) | deepplant_rebed_cpg
        attachReferenceFiles(deepplant_aggregate_chg.output, reference_map_ch) | deepplant_rebed_chg
        attachReferenceFiles(deepplant_aggregate_chh.output, reference_map_ch) | deepplant_rebed_chh

        deepplant_consolidate_cpg(deepplant_rebed_cpg.output)
        deepplant_consolidate_chg(deepplant_rebed_chg.output)
        deepplant_consolidate_chh(deepplant_rebed_chh.output)
}

