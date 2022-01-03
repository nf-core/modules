#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_AMPLICONCLIP } from '../../../../modules/samtools/ampliconclip/main.nf'

workflow test_samtools_ampliconclip_no_stats_no_rejects {

    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    save_cliprejects = false
    save_clipstats   = false

    SAMTOOLS_AMPLICONCLIP ( input, bed, save_cliprejects, save_clipstats )
}

workflow test_samtools_ampliconclip_no_stats_with_rejects {

    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    save_cliprejects = true
    save_clipstats   = false

    SAMTOOLS_AMPLICONCLIP ( input, bed, save_cliprejects, save_clipstats )
}

workflow test_samtools_ampliconclip_with_stats_with_rejects {

    input = [ 
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)
    save_cliprejects = true
    save_clipstats   = true

    SAMTOOLS_AMPLICONCLIP ( input, bed, save_cliprejects, save_clipstats )
}
