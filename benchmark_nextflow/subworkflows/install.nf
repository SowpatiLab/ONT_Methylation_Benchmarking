include {
    install_rerio;
    install_dorado ;
} from '../modules/dorado.nf'

include { setup_rockfish } from '../modules/rockfish.nf'
include { setup_f5c } from '../modules/f5c.nf'
include { setup_deepmod2 } from '../modules/deepmod2.nf'

workflow INSTALL_ALL_TOOLS {
    main:
        install_rerio()
        install_dorado(Channel.fromList(['dorado-0.9.1', 'dorado-1.1.1']))
        setup_rockfish()
        setup_f5c()
        setup_deepmod2()
}