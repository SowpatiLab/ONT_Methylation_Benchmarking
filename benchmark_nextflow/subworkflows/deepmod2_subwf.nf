include {
    setup_deepmod2 ;
    deepmod2_call;
    deepmod2_add_ref;
    consolidate_deepmod2
} from '../modules/deepmod2.nf'

// include {
//     setup_deepmod2 ;
//     deepmod2_call as deepmod2_call_transformer ;
//     deepmod2_add_ref as deepmod2_add_ref_transformer ;
//     consolidate_deepmod2 as consolidate_deepmod2_transformer
// } from '../modules/deepmod2.nf'
// include {
//     deepmod2_call as deepmod2_call_bilstm ;
//     deepmod2_add_ref as deepmod2_add_ref_bilstm ;
//     consolidate_deepmod2 as consolidate_deepmod2_bilstm
// } from '../modules/deepmod2.nf'

workflow DEEPMOD2_WORKFLOW {
    take: 
        samtools_cleansed_ch
        reference_map_ob
        model_type // 'transformer' or 'bilstm'
    main:
        setup_deepmod2()

        // remap pod5 files back to their bam files
        deepmod2_input = samtools_cleansed_ch
                .map{ bam, csi -> ["${bam.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", bam, csi ]}
                .combine(reference_map_ob, by: 0)
                .map{ it -> [ *it[1..-3],  file(it[-1])] }
                .combine(setup_deepmod2.output)
                .map{ it -> [ it[-1], *it[0..-2]] }

        deepMod2_modele_selector = deepmod2_input.combine(Channel.of(model_type))
        deepmod2_call(deepMod2_modele_selector)

        deemod2_calls_tuple = deepmod2_call.out.folder
            .combine(Channel.of(model_type))

        deepmod2_add_ref(deemod2_calls_tuple)
        consolidate_deepmod2(deepmod2_add_ref.output)
    // emit:
}
