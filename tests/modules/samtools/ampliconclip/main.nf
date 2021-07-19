#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_AMPLICONCLIP } from '../../../../modules/samtools/ampliconclip/main.nf' addParams([:])

workflow test_samtools_ampliconclip_no_stats_no_rejects {
    params.save_cliprejects = false
    params.save_clipstats   = false

    input = [ [ id:'test', single_end:false ],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    SAMTOOLS_AMPLICONCLIP ( input, bed )
}

workflow test_samtools_ampliconclip_no_stats_with_rejects {
    params.save_cliprejects = true
    params.save_clipstats   = false

    input = [ [ id:'test', single_end:false ],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    SAMTOOLS_AMPLICONCLIP ( input, bed )
}

workflow test_samtools_ampliconclip_with_stats_with_rejects {
    params.save_cliprejects = true
    params.save_clipstats   = true

    input = [ [ id:'test', single_end:false ],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]
    bed = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    SAMTOOLS_AMPLICONCLIP ( input, bed )
}
