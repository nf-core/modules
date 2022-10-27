#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { TIDDIT_COV as TIDDIT_COV_BED } from "$moduleDir/modules/nf-core/tiddit/cov/main.nf"
include { TIDDIT_COV as TIDDIT_COV_WIG } from "$moduleDir/modules/nf-core/tiddit/cov/main.nf"

workflow test_tiddit_cov_cram_bed {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    TIDDIT_COV_BED ( input, fasta )
}

workflow test_tiddit_cov_bam_bed {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    TIDDIT_COV_BED ( input, [] )
}

workflow test_tiddit_cov_cram_wig {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true) ]

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    TIDDIT_COV_WIG ( input, fasta )
}

workflow test_tiddit_cov_bam_wig {

    input = [ [ id:'test', single_end:false ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]

    TIDDIT_COV_WIG ( input, [] )
}
