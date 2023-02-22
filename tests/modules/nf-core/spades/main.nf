#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../../../../modules/nf-core/spades/main.nf'

workflow test_spades_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
              [],
              []
            ]
    SPADES ( input, [] , [] )
}

workflow test_spades_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
              [],
              []
            ]

    SPADES ( input, [] , [] )
}

workflow test_spades_illumina_nanopore {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
              [],
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]
            ]

    SPADES ( input, [] , [] )
}


workflow test_spades_yml {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
              [],
              [ file(params.test_data['sarscov2']['nanopore']['test_fastq_gz'], checkIfExists: true) ]

            ]
    yml = [
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/spades/spades_input_yml.yml', checkIfExists: true)
        ]
    SPADES ( input, yml, [] )
}

// that isnt perfect, because CCS reads should rather be used with -s instead of --pacbio
workflow test_spades_illumina_pacbio {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['test_2_fastq_gz'], checkIfExists: true) ],
              [ file(params.test_data['homo_sapiens']['pacbio']['ccs_fq_gz'], checkIfExists: true) ],
              []
            ]

    SPADES ( input, [] , [] )
}


