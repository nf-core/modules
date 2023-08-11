#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_DNASCOPE      } from '../../../../../modules/nf-core/sentieon/dnascope/main.nf'
include { SENTIEON_DNAMODELAPPLY } from '../../../../../modules/nf-core/sentieon/dnamodelapply/main.nf'

workflow test_dnamodelapply {

    fasta = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]).collect()
    fai   = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]).collect()

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    ml_model = file("https://s3.amazonaws.com/sentieon-release/other/SentieonDNAscopeModel1.0.model", checkIfExists: true)

    SENTIEON_DNASCOPE ( input, fasta, fai, [[:],[]], [[:],[]], [[:],[]], ml_model )

    ch_applyin = SENTIEON_DNASCOPE.out.vcf.join(SENTIEON_DNASCOPE.out.index)

    SENTIEON_DNAMODELAPPLY (ch_applyin, fasta, fai, ml_model)
}
