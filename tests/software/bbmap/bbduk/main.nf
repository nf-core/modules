#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BBMAP_BBDUK } from '../../../../software/bbmap/bbduk/main.nf' addParams( options: [:] )

workflow test_bbmap_bbduk_single_end {
    
    input = [ [ id:'test', single_end:true ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    contaminants = file('contaminants_dummy')
    use_contaminants = false

    BBMAP_BBDUK ( input, contaminants, use_contaminants )
}

workflow test_bbmap_bbduk_paired_end {
    
    input = [ [ id:'test', single_end:false ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                 file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    contaminants = file('contaminants_dummy')
    use_contaminants = false

    BBMAP_BBDUK ( input, contaminants, use_contaminants )
}

workflow test_bbmap_bbduk_se_ref {
    
    input = [ [ id:'test', single_end:true ], // meta map
              [  file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
                 
            ]
    contaminants = [file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ] // transciptome file - remove contaminants (*trim.fastq files empty)
    use_contaminants = true

    BBMAP_BBDUK ( input, contaminants, use_contaminants )
}


workflow test_bbmap_bbduk_pe_ref {
    
    input =  [  [ id:'test', single_end:false ], // meta map
                [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                  file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
             ]
    contaminants = [file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]
    use_contaminants = true

    BBMAP_BBDUK ( input, contaminants, use_contaminants )
}
