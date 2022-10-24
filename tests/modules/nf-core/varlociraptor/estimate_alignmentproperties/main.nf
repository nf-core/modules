#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATE_ALIGNMENTPROPERTIES } from '../../../../../modules/nf-core/varlociraptor/estimate_alignmentproperties/main.nf'

workflow test_varlociraptor_estimate_alignmentproperties {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    fasta = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai = [
        [ id:'sarscov2' ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_ESTIMATE_ALIGNMENTPROPERTIES ( input, fasta, fai )
}
