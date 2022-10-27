#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { BWA_INDEX                    } from "$moduleDir/modules/nf-core/bwa/index/main.nf"
include { TIDDIT_SV                    } from "$moduleDir/modules/nf-core/tiddit/sv/main.nf"
include { TIDDIT_SV as TIDDIT_SV_NOBWA } from "$moduleDir/modules/nf-core/tiddit/sv/main.nf"

workflow test_tiddit_sv_bam {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX( fasta )

    TIDDIT_SV ( input, fasta, BWA_INDEX.out.index)
}

workflow test_tiddit_sv_cram {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ],
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BWA_INDEX( fasta )

    TIDDIT_SV ( input, fasta, BWA_INDEX.out.index)
}

workflow test_tiddit_sv_nobwa_bam {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ],
        [ file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    TIDDIT_SV_NOBWA ( input, fasta, [])
}

workflow test_tiddit_sv_nobwa_cram {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ],
        [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true) ]
    ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    TIDDIT_SV_NOBWA ( input, fasta, [])
}
