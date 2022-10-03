#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ELPREP_FILTER } from '../../../../../modules/nf-core/elprep/filter/main.nf'

workflow test_elprep_filter {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    reference_elfasta = file(params.test_data['homo_sapiens']['genome']['genome_elfasta'], checkIfExists: true)
    known_sites_elsites = file(params.test_data['homo_sapiens']['genome']['dbsnp_146_hg38_elsites'], checkIfExists: true)
    target_regions_bed = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    ELPREP_FILTER ( input, true, true, [],  [], reference_elfasta, known_sites_elsites, target_regions_bed, [], [], true, true)
}
