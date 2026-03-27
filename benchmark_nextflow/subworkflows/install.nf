include {
    select_base_model;
    select_mod_model;
    fetch_ref ;
} from '../lib/utils.groovy'

include { setup_rockfish } from '../modules/rockfish.nf'
include { setup_f5c } from '../modules/f5c.nf'
include { setup_deepmod2 } from '../modules/deepmod2.nf'

workflow INSTALL_ALL_TOOLS {
    main:
        setup_rockfish()
        setup_f5c()
        setup_deepmod2()
}