#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATE_ALIGNMENTPROPERTIES } from '../../../../../modules/nf-core/varlociraptor/estimate_alignmentproperties/main.nf'
include { VARLOCIRAPTOR_PREPROCESS } from '../../../../../modules/nf-core/varlociraptor/preprocess/main.nf'

workflow test_varlociraptor_preprocess {

    bam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    bai = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    fasta = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    candidates = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_bcf'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_ESTIMATE_ALIGNMENTPROPERTIES ( bam, fasta, fai )
    VARLOCIRAPTOR_PREPROCESS ( bam, bai , VARLOCIRAPTOR_ESTIMATE_ALIGNMENTPROPERTIES.out.json, candidates, fasta, fai )
}
