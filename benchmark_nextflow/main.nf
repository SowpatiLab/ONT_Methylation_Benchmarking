include {
    select_base_model;
    select_mod_model;
} from './lib/utils.groovy'

include {
    install_rerio;
    download_rerio_models ;
    install_dorado ;
    download_dorado_base_model ;
    download_dorado_mod_models ;
    dorado_mod ;
    dorado_move;
} from './modules/dorado.nf'

include { samtools_move_cleanse ; samtools_mod_cleanse } from './modules/samtools.nf'
include { install_modkit; modkit_pileup ; modkit_add_ref ; standardise_dorado } from './modules/modkit.nf'
include {
    setup_rockfish ;
    rockfish_call ;
    rockfish_map_generate ;
    rockfish_intersect ;
    rockfish_aggregate ;
    rockfish_getfasta ;
    consolidate_rockfish
} from './modules/rockfish.nf'

include {
    setup_deepmod2 ;
    deepmod2_call as deepmod2_call_transformer ;
    deepmod2_add_ref as deepmod2_add_ref_transformer ;
    consolidate_deepmod2 as consolidate_deepmod2_transformer
} from './modules/deepmod2.nf'
include {
    deepmod2_call as deepmod2_call_bilstm ;
    deepmod2_add_ref as deepmod2_add_ref_bilstm ;
    consolidate_deepmod2 as consolidate_deepmod2_bilstm
} from './modules/deepmod2.nf'
include { bam_to_fastq ; pod5_to_blow5 ; f5c_idx_fastq ; f5c_call ; f5c_restrand } from './modules/f5c.nf'
include { setup_f5c ; f5c_aggregate ; f5c_add_fasta ; consolidate_f5c } from './modules/f5c.nf'
include {
    f5c_aggregate as f5c_aggregate_stranded ;
    f5c_add_fasta as f5c_add_fasta_stranded ;
    consolidate_f5c as consolidate_f5c_stranded
} from './modules/f5c.nf'

include { bam_fn_reorder ; deepbam_call ; deepplant_call } from './modules/deeptools.nf'
include {
    deetool_aggregate as deepbam_aggregate ;
    deeptool_rebed as deepbam_rebed ;
    deepbam_consolidate
} from './modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_cpg ; deeptool_rebed as deepplant_rebed_cpg ; deepplant_consolidate as deepplant_consolidate_cpg } from './modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_chg ; deeptool_rebed as deepplant_rebed_chg ; deepplant_consolidate as deepplant_consolidate_chg } from './modules/deeptools.nf'
include { deetool_aggregate as deepplant_aggregate_chh ; deeptool_rebed as deepplant_rebed_chh ; deepplant_consolidate as deepplant_consolidate_chh } from './modules/deeptools.nf'

include {
    nanostat;
    nanoplot;
    mosdepth;
} from './modules/qc.nf'


// includeConfig

workflow {

    reference_map = new groovy.yaml.YamlSlurper().parseText(file(params.refLookup).text)
    mod_list_filter = params.runControls.dorado_mod_models.findAll { e -> !e.value.isEmpty() }

    def run_qc = params.runControls.qc.tools.values().findAll{ it==true }.size() > 0 ? 0b01 : 0b00
    def others = params.runControls.other_tools.values().flatten().size() > 0 ? 0b10 : 0b00
    def run_selector =  run_qc + others
    // install dorado if not exist
    // check required models and determine which dorado versions to download

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
    install_versions = model_setup.map { dorver, m, v -> dorver }.unique()
    install_dorado(install_versions) // install dorado versions

    dorver_nam_aliases = params.toolConfig.dorado.dorado_version_aliases.collectEntries { key, value -> [(value): key] }

    install_dorado_versions = install_dorado.output
        .map { d -> [dorver_nam_aliases["$d.baseName"], d] }

    run_all      = model_setup
                    .combine(install_dorado_versions, by: 0)
                    
    run_mod  = run_all.filter { i -> i[1]==~/^[456]m[ACG]+$/ }
    run_move = run_all.filter { i -> !(i[1]==~/^[456]m[ACG]+$/) }
                    .map{ dorall, mod, ver, dorado -> 
                        def modelType = mod.matches('^[456]m[ACG]+$') ? mod : '5mC'
                        // def toolType = mod.matches('^[456]m[ACG]+$') ? 'dorado' : mod
                        
                        [ dorall, modelType, ver, dorado ]
                    }.unique()

    dorado_mod_run_channels = run_mod
        .combine(Channel.fromList(params.runControls.experiments))
        .combine(Channel.fromList(params.runControls.accuracy))
        .combine(Channel.of('5'))
        .unique()
        .map{ doralias, mod, ver, dorado, exp, acc, sr -> ["${sr}kHz_${acc}_${ver}_${mod}", exp ] }

    // setting up the input for base model download
    base_models = run_all
        .map{ dorall, mod, ver, dorado -> 
            def modelType = mod.matches('^[456]m[ACG]+$') ? mod : '5mC'            
            [ dorall, modelType, ver, dorado ]
        }
        .combine(Channel.fromList(params.runControls.experiments))
        .combine(Channel.fromList(params.runControls.accuracy))
        .combine(Channel.of('5'))
        .unique()
        .map{ doralias, mod, ver, dorado, exp, acc, sr -> [ select_base_model(sr, acc, ver, mod), dorado ] }
        .unique()

    download_dorado_base_model(base_models) // download base models
    
    if (mod_list_filter.size()>0){
        // setting up the input for mod model download
        mod_models = run_mod
            .combine(Channel.fromList(params.runControls.experiments))
            .combine(Channel.fromList(params.runControls.accuracy))
            .combine(Channel.of('5'))
            .unique()
            .map{ doralias, mod, ver, dorado, exp, acc, sr -> [ select_base_model(sr, acc, ver, mod), "${sr}kHz_${acc}_${ver}_${mod}", select_mod_model(sr, acc, ver, "_${mod}"), dorado ] }
            .unique()

        // filter out non existant models
        def black_list_models = params.runControls.dorado_model_blacklist
        mod_models = mod_models
            .filter{ it ->  
                if(black_list_models.contains(it[1].toString()) ){
                    log.warn "${it[0].toString()} not available, skipping."
                }
                !black_list_models.contains(it[1].toString()) }
            .branch{ it -> 
                rerio:  it[2].find('^res_*')
                dorado: it[2].find('^dna_*')
            }
        
        download_dorado_mod_models(mod_models.dorado)  // downalod dorado modification models
        
        def rerio_models = Channel.empty()

        mod_models.rerio.ifEmpty([])
        install_rerio()  // downalod rerio modification models
        rerio_models = mod_models.rerio
            .map{ baseAlias, modAlias, mod, dorado ->  [ baseAlias, modAlias, mod, dorado ]}
            .combine(install_rerio.output)

        rerio_models = download_rerio_models(rerio_models)
        
        all_mod_models =  rerio_models
            .mix(download_dorado_mod_models.output)

        filter_models = Channel.fromList(params.runControls.dorado_mod_models.values().flatten())
            .combine(Channel.fromList(params.runControls.accuracy))
            .map{ ver, acc -> [ select_base_model('5', acc, ver, '5mC'), 'dorado' ]}
            .unique()

        // remove all versions not listed in dorado_mod_models
        download_dorado_base_model_selected = download_dorado_base_model.output
            .combine(filter_models)
            .filter{ it -> it[0]==it[3]}
            .map{ it -> it[0..2]}

        dorado_mod_inputs = download_dorado_base_model_selected
            .combine(all_mod_models, by: 0)
            .map{ model_alias, dorado, base, key, mod -> [ key, dorado, base, mod ]}
            .combine(dorado_mod_run_channels, by: 0)
            .unique()
            .map{ key, dorver, base, mod, exp -> [exp, key, dorver, base, mod]}

        // mod basecall stage

        dorado_mod(dorado_mod_inputs)
        samtools_mod_cleanse(dorado_mod.output)

        install_modkit()
        mokit_pileup_in = samtools_mod_cleanse.out.bam
            .merge(samtools_mod_cleanse.out.csi)

        mokit_pileup_in = install_modkit.output
            .combine(mokit_pileup_in)
            

        modkit_pileup(mokit_pileup_in)
        modkit_add_ref(modkit_pileup.output)
        standardise_dorado(modkit_add_ref.output)
    }

    // todo: edge case, different model combinations across other tools
    /*
        OTHER TOOLS MODIFIED BASECALLING SECTION
    */
    // run_selector = 0
    if (run_selector > 0){

        def fastq = Channel.empty()
        
        filter_models = Channel.fromList(params.runControls.other_tools.values().flatten())
                .combine(Channel.fromList(params.runControls.accuracy))
                .map{ ver, acc -> [ select_base_model('5', acc, ver, '5mC'), 'dorado' ]}
                .unique()
                .ifEmpty([])
        
        // collect dorado models        
        download_dorado_base_model_move = download_dorado_base_model.output
                .combine(filter_models)
                .filter{ it -> it[0]==it[3]}
                .map{ it -> it[0..2]}
                
        // collect experiments to be run

        dorado_move_run_channels = run_move
            .unique()
            .combine(Channel.fromList(params.runControls.experiments))
            .combine(Channel.fromList(params.runControls.accuracy))
            .combine(Channel.of('5'))
            .unique()
            .map{ doralias, mod, ver, dorado, exp, acc, sr -> [ select_base_model(sr, acc, ver, mod), "${sr}kHz_${acc}_${ver}", exp ] }

        //
        dorado_move_inputs = download_dorado_base_model_move
            .map{ model_alias, dorado, base -> [ model_alias, dorado, base ]}
            .combine(dorado_move_run_channels, by: 0)
            .unique()
            .map{ model_alias, dorver, base, key, exp -> [ exp, key, dorver, base ]}
        
        dorado_move(dorado_move_inputs)
        samtools_move_cleanse(dorado_move.output)
        samtools_tuppled = samtools_move_cleanse.out.bam.merge( samtools_move_cleanse.out.csi )

        def run_f5c = ( params.runControls.other_tools.f5c.size() > 0 || params.runControls.other_tools.f5c_stranded.size() > 0)
        
        if (run_selector==1 || run_f5c ){
            bam_to_fastq(samtools_tuppled)
            fastq = bam_to_fastq.output
        }

        // /* ROCKFISH */
        if(params.runControls.other_tools.rockfish.size()>0){
            setup_rockfish()

            rockfish_input = setup_rockfish.output
                .combine(samtools_tuppled)

            rockfish_call(rockfish_input)

            rockfish_map_generate(samtools_tuppled.combine(rockfish_call.output))
            rockfish_intersect(rockfish_map_generate.output)
            rockfish_aggregate(rockfish_intersect.output)
            rockfish_getfasta(rockfish_aggregate.output)
            consolidate_rockfish(rockfish_getfasta.output)
        }

        /* DEEPOMOD2 */
        
        if (params.runControls.other_tools.deepmod2_transformer.size()>0 || params.runControls.other_tools.deepmod2_bilstm.size()>0){

            setup_deepmod2()

            deepmod2_input = setup_deepmod2.output.combine(samtools_tuppled)

            /* DEEPOMOD2 TRANSFORMER*/
            if (params.runControls.other_tools.deepmod2_transformer.size()>0){
                dm2_trans_channel = deepmod2_input.combine(Channel.of('transformer'))
                deepmod2_call_transformer(dm2_trans_channel)

                deemod2_calls_transformer_tuple = deepmod2_call_transformer.out.folder
                    .combine(Channel.of('transformer'))

                deepmod2_add_ref_transformer(deemod2_calls_transformer_tuple)
                consolidate_deepmod2_transformer(deepmod2_add_ref_transformer.output)
            }

            /* DEEPOMOD2 BiLSTM*/
            if (params.runControls.other_tools.deepmod2_bilstm.size()>0){
                dm2_blstm_channel = deepmod2_input.combine(Channel.of('BiLSTM'))
                deepmod2_call_bilstm(dm2_blstm_channel)

                deemod2_calls_bilstm_tuple = deepmod2_call_bilstm.out.folder
                    .combine(Channel.of('BiLSTM'))

                deepmod2_add_ref_bilstm(deemod2_calls_bilstm_tuple)
                consolidate_deepmod2_bilstm(deepmod2_add_ref_bilstm.output)
            }
        }

        // /* F5C */
        if( run_f5c ){            
            pod5_to_blow5(fastq)

            // installl f5c if not exist
            setup_f5c()

            f5c_idx_input  = setup_f5c.output.combine(pod5_to_blow5.output)
            f5c_idx_fastq(f5c_idx_input)
            f5c_call(f5c_idx_fastq.output)

            // /* F5C unstranded pipeline */
            if(params.runControls.other_tools.f5c.size()>0){
                f5c_call_ob  = setup_f5c.output.combine(f5c_call.output)
                f5c_aggregate(f5c_call_ob)
                f5c_add_fasta(f5c_aggregate.output)
                consolidate_f5c(f5c_add_fasta.output)
            }

            // // /* F5C restrand pipeline */
            if(params.runControls.other_tools.f5c_stranded.size()>0){
                f5c_restrand(f5c_call.output)
                f5c_call_ob  = setup_f5c.output.combine(f5c_restrand.output)
                f5c_aggregate_stranded(f5c_call_ob)
                f5c_add_fasta_stranded(f5c_aggregate_stranded.output)
                consolidate_f5c_stranded(f5c_add_fasta_stranded.output)
            }
        }
    
        /* DEEPBAM */
        if(params.runControls.other_tools.deepbam.size()>0 || params.runControls.other_tools.deepplant.size()>0 ){
            bam_fn_reorder(samtools_tuppled)
        }

        if(params.runControls.other_tools.deepbam.size()>0){
            deepbam_call(bam_fn_reorder.output)

            deepbam_aggregate(deepbam_call.output)
            deepbam_rebed(deepbam_aggregate.output)
            deepbam_consolidate(deepbam_rebed.output)
        }

        // /* DEEP_PLANT */
        if(params.runControls.other_tools.deepplant.size()>0){
            deepplant_call(bam_fn_reorder.output)
            
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

        if(run_qc){
            if(params.runControls.qc.tools.nanostat==true) { nanostat(samtools_tuppled) }
            if(params.runControls.qc.tools.nanoplot==true) { nanoplot(samtools_tuppled) }
            if(params.runControls.qc.tools.mosdepth==true) { mosdepth(samtools_tuppled) }
        }
    }
}
