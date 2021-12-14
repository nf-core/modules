#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBDUK } from '../../../../modules/bbmap/bbduk/main.nf'

workflow test_bbmap_bbduk_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    BBMAP_BBDUK ( input, [] )
}

workflow test_bbmap_bbduk_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    BBMAP_BBDUK ( input, [] )
}

workflow test_bbmap_bbduk_se_ref {

    input = [ [ id:'test', single_end:true ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    contaminants = [file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ] // transciptome file - remove contaminants (*trim.fastq files empty)

    BBMAP_BBDUK ( input, contaminants )
}

workflow test_bbmap_bbduk_pe_ref {

    input =  [  [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
             ]
    contaminants = [file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]

    BBMAP_BBDUK ( input, contaminants )
}
