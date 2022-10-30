#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SVTK_RDTEST2VCF } from '../../../../../modules/nf-core/svtk/rdtest2vcf/main.nf'
include { SVTK_VCF2BED      } from '../../../../../modules/nf-core/svtk/vcf2bed/main.nf'
include { MANTA_GERMLINE    } from '../../../../../modules/nf-core/manta/germline/main.nf'

workflow test_svtk_rdtest2vcf {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fasta_fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    MANTA_GERMLINE (
        input,
        fasta,
        fasta_fai
    )

    SVTK_VCF2BED (
        MANTA_GERMLINE.out.diploid_sv_vcf.combine(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by:0)
    )

    samples = Channel.of("test").collectFile(name: "samples.txt")

    SVTK_RDTEST2VCF (
        SVTK_VCF2BED.out.bed.combine(samples),
        fasta_fai
    )
}
