#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_SVCLUSTER } from '../../../../../modules/nf-core/gatk4/svcluster/main.nf'
include { MANTA_GERMLINE  } from '../../../../../modules/nf-core/manta/germline/main.nf'

workflow test_gatk4_svcluster {

    input = Channel.of([
        [ id:'normal' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ],
    [
        [ id:'tumor' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_cram_crai'], checkIfExists: true),
        [],
        []
    ])

    fasta     = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)])
    fasta_fai = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)])
    dict      = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)])

    MANTA_GERMLINE(
        input,
        fasta,
        fasta_fai
    )

    svcluster_input = MANTA_GERMLINE.out.diploid_sv_vcf.combine(
                        MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by: 0
                    ).map({ meta, vcf, tbi -> [ [id:'test'], vcf, tbi ]}).groupTuple()

    ploidy = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/svcluster/samples_ploidy.tsv", checkIfExists:true)

    GATK4_SVCLUSTER (
        svcluster_input,
        ploidy,
        fasta.map{ meta, fasta -> [fasta] },
        fasta_fai.map{ meta, fasta_fai -> [fasta_fai] },
        dict.map{ meta, dict -> [dict] }
    )
}
