#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX } from '../../../../../modules/nf-core/samtools/faidx/main.nf'

workflow test_samtools_faidx {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
              []
            ]

    SAMTOOLS_FAIDX ( input )
}

workflow test_samtools_faidx_bgzip {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true),
              []
            ]

    SAMTOOLS_FAIDX ( input )
}

workflow test_samtools_faidx_fasta {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
            ]

    SAMTOOLS_FAIDX ( input )
}

workflow test_samtools_faidx_stub_fasta {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
            ]

    SAMTOOLS_FAIDX ( input )
}

workflow test_samtools_faidx_stub_fai {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
              []
            ]

    SAMTOOLS_FAIDX ( input )
}
