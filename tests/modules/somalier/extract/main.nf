#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT } from '../../../../modules/somalier/extract/main.nf'

input_sample_ch = Channel
.fromPath( params.inputCsv )
.splitCsv( header:true )
.map { row -> tuple( [id: row.patient], file(row.bam), file(row.bai) ) }

fasta       = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta'], checkIfExists: true)

fasta_fai   = file(params.test_data['homo_sapiens']['genome']['genome_21_fasta_fai'], checkIfExists: true)

sites       = file(params.sites, checkIfExists: true)


workflow test_somalier_extract {

    SOMALIER_EXTRACT ( input_sample_ch, fasta, fasta_fai, sites )
}
