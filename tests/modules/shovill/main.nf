#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHOVILL                    } from '../../../modules/shovill/main.nf' addParams( options: [args: '--gsize 2800000 --kmers 31']          )
include { SHOVILL as SHOVILL_SKESA   } from '../../../modules/shovill/main.nf' addParams( options: [args: '--assembler skesa --gsize 2800000']   )
include { SHOVILL as SHOVILL_MEGAHIT } from '../../../modules/shovill/main.nf' addParams( options: [args: '--assembler megahit --gsize 2800000'] )
include { SHOVILL as SHOVILL_VELVET  } from '../../../modules/shovill/main.nf' addParams( options: [args: '--assembler velvet --gsize 2800000']  )

workflow test_shovill {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("https://github.com/nf-core/test-datasets/raw/bacass/ERR044595_1M_1.fastq.gz", checkIfExists: true),
                file("https://github.com/nf-core/test-datasets/raw/bacass/ERR044595_1M_2.fastq.gz", checkIfExists: true) ]
            ]

    SHOVILL ( input )
}

workflow test_shovill_megahit {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    SHOVILL_MEGAHIT ( input )
}

workflow test_shovill_skesa {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    SHOVILL_SKESA ( input )
}

workflow test_shovill_velvet {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    SHOVILL_VELVET ( input )
}
