#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWA_INDEX } from '../../../../modules/bwa/index/main.nf'
include { TIDDIT_SV } from '../../../../modules/tiddit/sv/main.nf'

workflow test_tiddit_sv {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)

    BWA_INDEX( fasta )

    TIDDIT_SV ( input, fasta, fai , BWA_INDEX.out.index)
}

workflow test_tiddit_sv_no_ref {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    ]

    TIDDIT_SV ( input, [], [] )
}

workflow test_tiddit_sv_cram {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai   = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    TIDDIT_SV ( input, fasta, fai )
}
