#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_VIEW } from '../../../../../modules/nf-core/samtools/view/main.nf'

workflow test_samtools_view {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
                []
            ]

    SAMTOOLS_VIEW ( input, [[],[]], [] )
}

workflow test_samtools_view_cram {
    input = [ [ id: 'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_VIEW ( input, fasta, [] )
}

workflow test_samtools_view_convert {
    input = [ [ id: 'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                []
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_VIEW ( input, fasta, [] )
}

workflow test_samtools_view_index {
    input = [ [ id: 'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                []
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    SAMTOOLS_VIEW ( input, fasta, [] )
}

workflow test_samtools_view_filter {
    input = [ [ id: 'test' ], // meta map
                file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
                []
            ]
    fasta = [ [ id:'genome' ],
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
            ]

    qname = Channel.of("testN:2817", "testN:2814").collectFile(name: "readnames.list", newLine: true)

    SAMTOOLS_VIEW ( input, fasta, qname )
}

workflow test_samtools_view_stubs {
    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
                []
            ]

    SAMTOOLS_VIEW ( input, [[],[]], [] )
}
