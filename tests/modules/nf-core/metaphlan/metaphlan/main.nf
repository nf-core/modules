#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR               } from '../../../../../modules/nf-core/untar/main.nf'
include { METAPHLAN_METAPHLAN } from '../../../../../modules/nf-core/metaphlan/metaphlan/main.nf'
include { SAMTOOLS_VIEW       } from '../../../../../modules/nf-core/samtools/view/main.nf'

workflow test_metaphlan_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan4_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN_METAPHLAN ( input, UNTAR.out.untar.map{ it[1] } )
}

workflow test_metaphlan_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan4_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN_METAPHLAN ( input, UNTAR.out.untar.map{ it[1] } )
}

workflow test_metaphlan_fasta {

    input = [ [ id:'test', single_end:true], // meta map
              [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan4_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    METAPHLAN_METAPHLAN ( input, UNTAR.out.untar.map{ it[1] } )
}

workflow test_metaphlan_sam {

    input = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
                []
            ]

    db    = [ [], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/metaphlan4_database.tar.gz', checkIfExists: true) ]

    UNTAR ( db )
    SAMTOOLS_VIEW ( input, [[],[]], [])
    METAPHLAN_METAPHLAN ( SAMTOOLS_VIEW.out.sam, UNTAR.out.untar.map{ it[1] } )
}

