nextflow.enable.dsl = 2

include {
    TRIMMOMATIC as TRIMMOMATIC_SE
    TRIMMOMATIC as TRIMMOMATIC_PE
    TRIMMOMATIC
} from '../../../../modules/nf-core/trimmomatic/main.nf'

//
// Test with single-end data
//
workflow test_trimmomatic_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    TRIMMOMATIC_SE ( input )
}

//
// Test with paired-end data
//
workflow test_trimmomatic_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    TRIMMOMATIC_PE ( input )
}

//
// Failing test with no adaptor
//
workflow test_trimmomatic_no_adaptor {
    input = [ [ id:'test', single_end:false ], // meta map
             [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    TRIMMOMATIC ( input )
}
