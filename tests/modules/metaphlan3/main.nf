#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR          } from '../../../modules/untar/main.nf'
include { SAMTOOLS_VIEW  } from '../../../modules/samtools/view/main.nf'
include { METAPHLAN3     } from '../../../modules/metaphlan3/main.nf'

workflow test_metaphlan3_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN3 ( input, UNTAR.out.untar.map{ it[1] } )
}

workflow test_metaphlan3_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN3 ( input, UNTAR.out.untar.map{ it[1] } )
}

workflow test_metaphlan3_sam {

    input = [ [ id:'test'], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    SAMTOOLS_VIEW ( input, [] )
    METAPHLAN3 ( SAMTOOLS_VIEW.out.bam, UNTAR.out.untar.map{ it[1] } )
}

workflow test_metaphlan3_fasta {

    input = [ [ id:'test', single_end:true], // meta map
              [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN3 ( input, UNTAR.out.untar.map{ it[1] } )
}
