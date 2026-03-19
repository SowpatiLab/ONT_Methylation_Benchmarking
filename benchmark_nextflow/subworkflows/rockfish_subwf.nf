include {
    setup_rockfish ;
    download_rockfish_model ;
    rockfish_call ;
    rockfish_map_generate ;
    rockfish_intersect ;
    rockfish_aggregate ;
    rockfish_getfasta ;
    consolidate_rockfish
} from '../modules/rockfish.nf'

include { attachReferenceFiles } from '../lib/utils.groovy'

workflow ROCKFISH_SUBWORKFLOW {
    take: 
        reference_map_ch
        samtools_cleansed_ch
        reference_map_ob
    main:   
        rockfish_exe = Channel.empty()
        if(workflow.containerEngine!=null){
            rockfish_exe = Channel.of('rockfish')
        } else {
            setup_rockfish()
            rockfish_exe = setup_rockfish.output
        }
        download_rockfish_model(rockfish_exe)
        rockfish_model = download_rockfish_model.output

        rockfish_exe
            .combine(samtools_cleansed_ch)
            .map{ it -> [ file("${it[2]}").baseName.toString().split('_5kHz')[0] + '_5kHz' , *it ]}
            .join(reference_map_ob, by: 0)
            .combine(rockfish_model)
            .map{ it -> [ it[4], it[1], it[6], *it[2..3]] } |
            rockfish_call

        samtools_cleansed_ch
            .combine(rockfish_call.output)
            .map{ bam, csi, rfOut -> ["${bam.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", bam, csi, rfOut ] }
            .combine( reference_map_ob, by: 0 )
            .map{ it -> [ *it[1..3], it[-1]] } |
            rockfish_map_generate |
            rockfish_intersect |
            rockfish_aggregate

            attachReferenceFiles(rockfish_aggregate.output, reference_map_ch) | 
            rockfish_getfasta |
            consolidate_rockfish
}