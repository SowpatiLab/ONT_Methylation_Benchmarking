include { 
    DORADO_SETUP; DORADO_BASECALL_MODIFIED; 
    DORADO_BASECALL_MOVETABLE 
} from './subworkflows/dorado_basecall_subwf.nf'

include { ROCKFISH_SUBWORKFLOW } from './subworkflows/rockfish_subwf.nf'

include {
    DEEPMOD2_WORKFLOW as DEEPMOD2_WORKFLOW_TRANSFORMER;
    DEEPMOD2_WORKFLOW as DEEPMOD2_WORKFLOW_BILSTM;
} from './subworkflows/deepmod2_subwf.nf'

include { 
    F5C_GENERAL_WORKFLOW;
    F5C_UNSTRANDED_WORKFLOW;
    F5C_RESTRANDED_WORKFLOW;
} from './subworkflows/f5c_subwf.nf' 

include { DEEPTOOLS_GENERAL_SUBWORKFLOW; DEEPBAM_SUBWORKFLOW; DEEPPLANT_SUBWORKFLOW } from './subworkflows/deeptools_subwf.nf'

include {
    nanostat;
    nanoplot;
    mosdepth;
    nanoq;
} from './modules/qc.nf'


// includeConfig

workflow {

    reference_map = new groovy.yaml.YamlSlurper().parseText(file(params.refLookup).text)
    def reference_map_ob = Channel.fromList(params.runControls.experiments)
                .combine(Channel.of(params.pod5dir))
                .map{ p5, dir -> [ 
                    "${p5}_5kHz", 
                    file("${dir}/${p5}_5kHz"), 
                    file(reference_map['references'][p5.split('_')[0]])
                ]}

    mod_list_filter = params.runControls.dorado_mod_models.findAll { e -> !e.value.isEmpty() }

    def run_qc = params.runControls.qc.tools.values().findAll{ it==true }.size() > 0 ? 0b01 : 0b00 // return 1 if qc tools are activated
    def others = params.runControls.other_tools.values().flatten().size() > 0 ? 0b10 : 0b00        // return 2 f other tools are activated
    def run_selector =  run_qc + others

    models_lookup_mods = params.runControls.dorado_mod_models.collectMany { k, vs -> vs.collect { v -> [k, v] } } // fetch mod models to run
    models_lookup_move = params.runControls.other_tools.collectMany { k, vs -> vs.collect { v -> [k, v] } }       // fetch move models to run
    
    // if qc is activate d set default model to v5r3
    models_lookup_move = (run_qc==1 && models_lookup_move.size()==0) ? [['qc', params.runControls.qc.default_model]] : models_lookup_move
    model_lookup_all = models_lookup_mods + models_lookup_move  // combine both

    // create a model lookup for later
    model_setup = Channel
        .from(model_lookup_all)
        .map { m, v -> [params.toolConfig.dorado.dorado_version_map["5kHz_${v.split('r')[0]}"], m, v] }
        .unique()

    // isolate unique dorado versions that need to be installed
    experiment_list = Channel.fromList(params.runControls.experiments)
    accuracy_list   = Channel.fromList(params.runControls.accuracy)
    install_versions = model_setup.map { dorver, m, v -> dorver }.unique()
    dorver_nam_aliases = params.toolConfig.dorado.dorado_version_aliases.collectEntries { key, value -> [(value): key] }

    /* DORADO SETUP */
    // // install dorado if not exist
    // // check required models and determine which dorado versions to download
    DORADO_SETUP( experiment_list, accuracy_list, model_setup, install_versions, dorver_nam_aliases )
    
    run_mod = DORADO_SETUP.out.run_mod
    run_move = DORADO_SETUP.out.run_move
    dorado_mod_run_channels = DORADO_SETUP.out.dorado_mod_run_channels
    download_dorado_base_model = DORADO_SETUP.out.download_dorado_base_model

    // RUN DORADO BASECALLER for MODIFIED BASES
    if(models_lookup_mods.size() > 0){
        def basecall_model_list_mod = Channel.fromList(params.runControls.dorado_mod_models.values().flatten())
        
        DORADO_BASECALL_MODIFIED(
            run_mod,
            experiment_list,
            accuracy_list,
            basecall_model_list_mod,
            params.runControls.black_list_models,
            dorado_mod_run_channels,
            download_dorado_base_model
        )
    }

    if(run_selector>0){
        
        def basecall_model_list_base = Channel.fromList(params.runControls.other_tools.values().flatten())
            .unique()
            .ifEmpty { params.runControls.qc.default_model }

        def run_f5c     = params.runControls.other_tools.f5c.size() > 0
        def run_f5c_stranded =  params.runControls.other_tools.f5c_stranded.size()  > 0
        def run_rockfish = params.runControls.other_tools.rockfish.size()>0
        def run_deepMod2_transformer = params.runControls.other_tools.deepmod2_transformer.size()>0
        def run_deepMod2_bilstm = params.runControls.other_tools.deepmod2_bilstm.size()>0
        def run_deepBAM   = params.runControls.other_tools.deepbam.size() > 0
        def run_deepPlant = params.runControls.other_tools.deepplant.size() > 0

        def generate_fastq = (run_selector==1 || run_f5c )
        
        experiment_list.view()
        basecall_model_list_base.view()
        
        DORADO_BASECALL_MOVETABLE(
            run_move,
            experiment_list,
            accuracy_list,
            basecall_model_list_base,
            download_dorado_base_model,
            generate_fastq
        )

        samtools_cleansed = DORADO_BASECALL_MOVETABLE.out.samtools_cleansed
        fastq = DORADO_BASECALL_MOVETABLE.out.fastq

        /* ROCKFISH */
        if(run_rockfish){
            ROCKFISH_SUBWORKFLOW( samtools_cleansed, reference_map_ob )
        }

        /* DEEPOMOD2 */
        if( run_deepMod2_transformer || run_deepMod2_bilstm ){
            if(run_deepMod2_transformer) {
                DEEPMOD2_WORKFLOW_TRANSFORMER (
                    samtools_cleansed,
                    reference_map_ob,
                    'transformer'
                )
            }
            if(run_deepMod2_bilstm) {
                DEEPMOD2_WORKFLOW_BILSTM (
                    samtools_cleansed,
                    reference_map_ob,
                    'BiLSTM'
                )
            }
        }
        /* F5C */
        if(run_f5c || run_f5c_stranded){
            F5C_GENERAL_WORKFLOW(fastq, reference_map_ob)
            def f5c_setup = F5C_GENERAL_WORKFLOW.out.setup_f5c
            def f5c_call  = F5C_GENERAL_WORKFLOW.out.f5c_call

            if(run_f5c) {
                F5C_UNSTRANDED_WORKFLOW(f5c_setup, f5c_call)
            }
            if(run_f5c_stranded){
                F5C_RESTRANDED_WORKFLOW(f5c_setup, f5c_call)
            }
        }
        /* DeepTools */
        if(run_deepBAM || run_deepPlant){
            def bam_fn_reorder = DEEPTOOLS_GENERAL_SUBWORKFLOW( samtools_cleansed )

            /* DEEPBAM */
            if(run_deepBAM){
                def deepbam_model_channel = Channel.of(params.toolConfig.deepbam.model)

                DEEPBAM_SUBWORKFLOW(
                    bam_fn_reorder,
                    reference_map_ob,
                    deepbam_model_channel
                )
            }
            /* DEEP_PLANT */
            if(run_deepPlant){
                def deepplant_model_channel = Channel.of(params.toolConfig.deepplant.model)

                DEEPPLANT_SUBWORKFLOW(
                    bam_fn_reorder,
                    reference_map_ob,
                    deepplant_model_channel
                )
            }
        }
        
        if(run_qc>0){
            if(params.runControls.qc.tools.nanostat==true) { nanostat(samtools_cleansed) }
            if(params.runControls.qc.tools.nanoplot==true) { nanoplot(samtools_cleansed) }
            if(params.runControls.qc.tools.mosdepth==true) { mosdepth(samtools_cleansed) }
            if(params.runControls.qc.tools.nanoq==true) { nanoq(fastq) }
        }
    }
}
