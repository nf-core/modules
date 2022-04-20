#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXPANSIONHUNTER } from '../../../modules/expansionhunter/main.nf'
include { STRANGER } from '../../../modules/stranger/main.nf'


input = [ [ id:'test', gender:'male' ], // meta map
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
            file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        ]
fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
variant_catalog = file(params.test_data['homo_sapiens']['genome']['repeat_expansions'], checkIfExists: true)


workflow test_stranger {
    EXPANSIONHUNTER ( input, fasta, variant_catalog )
    STRANGER ( EXPANSIONHUNTER.out.vcf, variant_catalog )
}

workflow test_stranger_without_optional_variant_catalog {
    EXPANSIONHUNTER ( input, fasta, variant_catalog )
    STRANGER ( EXPANSIONHUNTER.out.vcf, [] )
}
