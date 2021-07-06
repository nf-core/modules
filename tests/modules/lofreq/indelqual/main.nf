#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { LOFREQ_INDELQUAL } from '../../../../software/lofreq/indelqual/main.nf' addParams( options: [ 'args': '--dindel', 'suffix':'.indelqual'] )

workflow test_lofreq_indelqual {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    LOFREQ_INDELQUAL ( input, fasta )
}
