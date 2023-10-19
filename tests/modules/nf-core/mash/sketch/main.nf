#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MASH_SKETCH } from '../../../../../modules/nf-core/mash/sketch/main.nf'

workflow test_mash_sketch {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
              ]
            ]

    MASH_SKETCH ( input, false, false )
}

workflow test_mash_sketch_reads {

    input = [ [ id:'test', single_end:false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
              ]
            ]

    MASH_SKETCH ( input, true, false )
}

workflow test_mash_sketch_multi_fasta {

    input = [ [ id:'test', single_end:false], // meta map
              [ file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)
              ]
            ]

    MASH_SKETCH ( input, false, true )
}
