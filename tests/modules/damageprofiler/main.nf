#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DAMAGEPROFILER } from '../../../modules/damageprofiler/main.nf' addParams( options: [:] )

workflow test_damageprofiler {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ] ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    DAMAGEPROFILER ( input, fasta, fai )
}
