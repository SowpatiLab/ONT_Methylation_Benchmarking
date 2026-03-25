include {
    install_rerio;
    download_rerio_models ;
    install_dorado ;
    download_dorado_base_model ;
    download_dorado_mod_models ;
    dorado_mod ;
    dorado_move;
} from '../modules/dorado.nf'

include {
    select_base_model;
    select_mod_model;
    fetch_ref ;
} from '../lib/utils.groovy'

include { samtools_move_cleanse ; samtools_mod_cleanse } from '../modules/samtools.nf'
include { install_modkit; modkit_pileup ; modkit_add_ref ; standardise_dorado } from '../modules/modkit.nf'
include { bam_to_fastq } from '../modules/f5c.nf'

include { attachReferenceFiles } from '../lib/utils.groovy'

workflow DORADO_SETUP {
    take:
        experiment_list_ch
        accuracy_list_ch
        model_setup_ch
        install_versionso_ch
        dorver_nam_aliases_ch

    main:
        if(workflow.containerEngine!=null){
            // if containerised
            install_dorado_versions = install_versionso_ch
                .map { d -> [d, d] }
        } else {
            // if not installed, install dorado
            install_dorado(install_versionso_ch)
            install_dorado_versions = install_dorado.output
                .map { d -> [dorver_nam_aliases_ch["$d.baseName"], d] }
        }

        run_all  = model_setup_ch
                        .combine(install_dorado_versions, by: 0)
                        
        run_mod  = run_all.filter { i -> i[1]==~/^[456]m[ACG]+$/ }
        run_move = run_all.filter { i -> !(i[1]==~/^[456]m[ACG]+$/) }
                        .map{ dorall, mod, ver, dorado -> 
                            def modelType = mod.matches('^[456]m[ACG]+$') ? mod : '5mC'
                            [ dorall, modelType, ver, dorado ]
                        }.unique()

        dorado_mod_run_channels = run_mod
            .combine(experiment_list_ch)
            .combine(accuracy_list_ch)
            .combine(Channel.of('5'))
            .unique()
            .map{ doralias, mod, ver, dorado, exp, acc, sr -> ["${sr}kHz_${acc}_${ver}_${mod}", exp ] }

        // setting up the input for base model download
        base_models = run_all
            .map{ dorall, mod, ver, dorado -> 
                def modelType = mod.matches('^[456]m[ACG]+$') ? mod : '5mC'            
                [ dorall, modelType, ver, dorado ]
            }
            .combine(experiment_list_ch)
            .combine(accuracy_list_ch)
            .combine(Channel.of('5'))
            .unique()
            .map{ doralias, mod, ver, dorado, exp, acc, sr -> [ select_base_model(sr, acc, ver, mod), dorado ] }
            .unique()

        download_dorado_base_model(base_models) // download base models

    emit:
        run_mod = run_mod
        run_move = run_move
        dorado_mod_run_channels = dorado_mod_run_channels
        download_dorado_base_model = download_dorado_base_model.output

}

workflow DORADO_BASECALL_MODIFIED {
    take:
        reference_map_ch
        run_mod_ch
        experiment_list_ch
        accuracy_list_ch
        full_model_list_ch
        black_list_models_ch
        dorado_mod_run_ch
        downloaded_base_models_ch
    main:
        // setting up the input for mod model download
        def IS_CONTAINERISED = workflow.containerEngine!=null

        mod_models = run_mod_ch
            .combine(experiment_list_ch)
            .combine(accuracy_list_ch)
            .combine(Channel.of('5'))
            .unique()
            .map{ doralias, mod, ver, dorado, exp, acc, sr -> [ select_base_model(sr, acc, ver, mod), "${sr}kHz_${acc}_${ver}_${mod}", select_mod_model(sr, acc, ver, "_${mod}"), dorado ] }
            .unique()

        // filter out non existant models
        mod_models = mod_models
            .filter{ it ->  
                if(black_list_models_ch.contains(it[1].toString()) ){
                    log.warn "${it[1].toString()} not available, skipping."
                }
                !black_list_models_ch.contains(it[1].toString()) }
            .branch{ it -> 
                rerio:  it[2].find('^res_*')
                dorado: it[2].find('^dna_*')
            }
        
        download_dorado_mod_models(mod_models.dorado)  // downalod dorado modification models
        
        def rerio_models = Channel.empty()

        mod_models.rerio.ifEmpty([])
        rerio_exec = Channel.empty()

        if(workflow.containerEngine!=null){
            rerio_exec = Channel.of('rerio')
        } else {
            install_rerio()
            rerio_exec = install_rerio.output  // downalod rerio modification models
        }

        rerio_models = mod_models.rerio
            .map{ baseAlias, modAlias, mod, dorado ->  [ baseAlias, modAlias, mod, dorado ]}
            .combine(rerio_exec)

        rerio_models = download_rerio_models(rerio_models)
        
        all_mod_models =  rerio_models
            .mix(download_dorado_mod_models.output)
        
        // accuracy_list_ch.view()

        filter_models = full_model_list_ch
            .combine(accuracy_list_ch)
            .map{ ver, acc -> [ select_base_model('5', acc, ver, '5mC'), 'dorado' ]}
            .unique()

        // remove all versions not listed in dorado_mod_models
        download_dorado_base_model_selected = downloaded_base_models_ch
            .combine(filter_models)
            .filter{ it -> it[0]==it[3]}
            .map{ it -> it[0..2]}

        dorado_mod_inputs = download_dorado_base_model_selected
            .combine(all_mod_models, by: 0)
            .map{ model_alias, dorado, base, key, mod -> [ key, dorado, base, mod ]}
            .combine(dorado_mod_run_ch, by: 0)
            .unique()
            .map{ key, dorver, base, mod, exp -> [exp, key, dorver, base, mod, file("${params.pod5dir}/${exp}_5kHz"),  file(fetch_ref(exp))]}

        // mod basecall stage
        dorado_mod(dorado_mod_inputs)
        samtools_mod_cleanse(dorado_mod.output)
        
        def modkit_exec = Channel.empty()
        if(IS_CONTAINERISED){
            modkit_exec = Channel.of('modkit')
        }else{
            install_modkit()
            modkit_exec = install_modkit.output
        }
        
        mokit_pileup_in = modkit_exec
            .combine(samtools_mod_cleanse.output)

        modkit_pileup(mokit_pileup_in)
        
        attachReferenceFiles(modkit_pileup.output, reference_map_ch) | 
            modkit_add_ref
            standardise_dorado(modkit_add_ref.output)
}


workflow DORADO_BASECALL_MOVETABLE {
    take:
        run_move_ch
        experiment_list_ch
        accuracy_list_ch
        basecall_model_list_base
        downloaded_base_models_ch
        generate_fastq_flag
    main:
        def fastq_ch = Channel.empty()

        filter_models = basecall_model_list_base
            .combine(accuracy_list_ch)
            .map{ ver, acc -> [ select_base_model('5', acc, ver, '5mC'), 'dorado' ]}
            .unique()
            .ifEmpty([])

        // collect dorado models        
        download_dorado_base_model_move = downloaded_base_models_ch
                .combine(filter_models)
                .filter{ it -> it[0]==it[3]}
                .map{ it -> it[0..2] }
                // .map{ it -> [*it[0..2], file("${params.pod5dir}/${exp}_5kHz"),  file(fetch_ref(exp))]}
        
        // collect experiments to be run
        dorado_move_run_channels = run_move_ch
            .unique()
            .combine(experiment_list_ch)
            .combine(accuracy_list_ch)
            .combine(Channel.of('5'))
            .unique()
            .map{ doralias, mod, ver, dorado, exp, acc, sr -> [ select_base_model(sr, acc, ver, mod), "${sr}kHz_${acc}_${ver}", exp ] }
        
        dorado_move_inputs = download_dorado_base_model_move
            .map{ model_alias, dorado, base -> [ model_alias, dorado, base ]}
            .combine(dorado_move_run_channels, by: 0)
            .unique()
            .map{ model_alias, dorver, base, key, exp -> [ exp, key, dorver, base, file("${params.pod5dir}/${exp}_5kHz"),  file(fetch_ref(exp)) ]}
        
        dorado_move(dorado_move_inputs)
        samtools_move_cleanse(dorado_move.output)

        if(generate_fastq_flag){
            bam_to_fastq(samtools_move_cleanse.output)
            fastq_ch = bam_to_fastq.output
        }

    emit:
        dorado_move = dorado_move.output
        samtools_cleansed = samtools_move_cleanse.output
        fastq = fastq_ch
}