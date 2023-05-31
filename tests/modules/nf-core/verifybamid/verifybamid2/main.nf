#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VERIFYBAMID_VERIFYBAMID2 as VERIFYBAMID_VERIFYBAMID2_SVD } from '../../../../../modules/nf-core/verifybamid/verifybamid2/main.nf'
include { VERIFYBAMID_VERIFYBAMID2 as VERIFYBAMID_VERIFYBAMID2_REFVCF } from '../../../../../modules/nf-core/verifybamid/verifybamid2/main.nf'
include { VERIFYBAMID_VERIFYBAMID2 as VERIFYBAMID_VERIFYBAMID2_PANEL } from '../../../../../modules/nf-core/verifybamid/verifybamid2/main.nf'

workflow test_verifybamid_verifybamid2_svd {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    svd_prefix = [ file(params.niche_test['test_genome']['test_genome_ud'], checkIfExists: true),
                   file(params.niche_test['test_genome']['test_genome_mu'], checkIfExists: true),
                   file(params.niche_test['test_genome']['test_genome_bed'], checkIfExists: true)
                 ]

    references = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    VERIFYBAMID_VERIFYBAMID2_SVD ( input, svd_prefix, [], references )
}

workflow test_verifybamid_verifybamid2_panel {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    svd_prefix = [ file(params.niche_test['test_genome']['test_genome_ud'], checkIfExists: true),
                   file(params.niche_test['test_genome']['test_genome_mu'], checkIfExists: true),
                   file(params.niche_test['test_genome']['test_genome_bed'], checkIfExists: true)
                 ]

    references = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    VERIFYBAMID_VERIFYBAMID2_PANEL ( input, svd_prefix, [], references )
}

workflow test_verifybamid_verifybamid2_refvcf {
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
            ]

    svd_prefix = [ file(params.niche_test['test_genome']['test_genome_ud'], checkIfExists: true),
                   file(params.niche_test['test_genome']['test_genome_mu'], checkIfExists: true),
                   file(params.niche_test['test_genome']['test_genome_bed'], checkIfExists: true)
                 ]

    refvcf = file(params.test_data['homo_sapiens']['illumina']['test_genome_vcf'], checkIfExists: true)

    references = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    VERIFYBAMID_VERIFYBAMID2_REFVCF ( input, [[],[],[]], refvcf, references )
}
