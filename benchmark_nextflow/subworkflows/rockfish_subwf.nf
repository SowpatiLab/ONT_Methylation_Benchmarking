include {
    setup_rockfish ;
    rockfish_call ;
    rockfish_map_generate ;
    rockfish_intersect ;
    rockfish_aggregate ;
    rockfish_getfasta ;
    consolidate_rockfish
} from '../modules/rockfish.nf'

workflow ROCKFISH_SUBWORKFLOW {
    take: 
        samtools_cleansed_ch
        reference_map_ob
    main:   
        setup_rockfish()

        setup_rockfish.output
            .combine(samtools_cleansed_ch)
            .map{ it -> [ file("${it[2]}").baseName.toString().split('_5kHz')[0] + '_5kHz' , *it ]}
            .join(reference_map_ob, by: 0)
            .map{ it -> [ it[5], *it[1..4]]} |
            rockfish_call

        samtools_cleansed_ch
            .combine(rockfish_call.output)
            .map{ bam, csi, rfOut -> ["${bam.toString().split('/')[-1].split('_5kHz')[0]}_5kHz", bam, csi, rfOut ] }
            .combine( reference_map_ob, by: 0 )
            .map{ it -> [ *it[1..3], it[-1]] } |
            rockfish_map_generate |
            rockfish_intersect |
            rockfish_aggregate |
            rockfish_getfasta |
            consolidate_rockfish
}