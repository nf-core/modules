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
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [] // no intervals
    ]

    pcr_indel_model = 'CONSERVATIVE'
    ml_model        = Channel.of([[:], file("https://s3.amazonaws.com/sentieon-release/other/SentieonDNAscopeModel1.1.model", checkIfExists: true)])
    emit_vcf        = 'variant'
    emit_gvcf       = false

    SENTIEON_DNASCOPE (
        input,
        fasta,
        fai,
        [[:], []],
        [[:], []],
        ml_model,
        pcr_indel_model, // pcr_indel_model
        emit_vcf,  // emit_vcf
        emit_gvcf) // emit_gvcf

    ch_applyin = SENTIEON_DNASCOPE.out.vcf.join(SENTIEON_DNASCOPE.out.vcf_tbi)

    SENTIEON_DNAMODELAPPLY (ch_applyin, fasta, fai, ml_model)
}
