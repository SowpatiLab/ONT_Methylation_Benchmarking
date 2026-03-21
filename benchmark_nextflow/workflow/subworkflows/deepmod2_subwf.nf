include {
    setup_deepmod2 ;
    deepmod2_call;
    deepmod2_add_ref;
    consolidate_deepmod2
} from '../modules/deepmod2.nf'

workflow DEEPMOD2_WORKFLOW {
    take: 
        reference_map_ch
        samtools_cleansed_ch
        reference_map_ob
        model_type // 'transformer' or 'bilstm'
    main:

        // deepmod2_exe = Channel.empty()
        setup_deepmod2()
        deepmod2_exe = setup_deepmod2.output
        // if(workflow.containerEngine!=null){
        //     deepmod2_exe = Channel.of('deepmod2')
        // } else {
        //     setup_deepmod2()
        //     deepmod2_exe = setup_deepmod2.output
        // }

        // remap pod5 files back to their bam files
        deepmod2_input = samtools_cleansed_ch
            .map{ bam, csi -> ["${bam.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", bam, csi ]}
            .combine(reference_map_ob, by: 0)
            .map{ it -> [ *it[1..-3],  file(it[-1])] }
            .combine(deepmod2_exe)
            .map{ it -> [ it[-1], *it[0..-2]] }

        deepMod2_modele_selector = deepmod2_input.combine(Channel.of(model_type))
        deepmod2_call(deepMod2_modele_selector)

        deepmod2_call.out.folder
            .combine(Channel.of(model_type))
            .map{ bed, model -> [ 
                bed, model, reference_map_ch['references'][bed.baseName.split('_')[0]],
            ]}
            .map{ bed, model, ref -> [
                bed,
                model,
                file(ref),
                file(ref.replaceFirst(/\.fa*$/, '.fa.fai')),
                file(ref.replaceFirst(/\.fa*$/, '.genome'))
            ]} |
            deepmod2_add_ref |
            consolidate_deepmod2
}
