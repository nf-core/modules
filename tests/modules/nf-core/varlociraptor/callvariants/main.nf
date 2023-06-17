#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../../../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main.nf'
include { VARLOCIRAPTOR_PREPROCESS                  } from '../../../../../modules/nf-core/varlociraptor/preprocess/main.nf'
include { VARLOCIRAPTOR_CALLVARIANTS                } from '../../../../../modules/nf-core/varlociraptor/callvariants/main.nf'

workflow test_varlociraptor_callvariants {

    bam = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fai= [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES( bam, fasta, fai)

    input = Channel.of([
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true),
    ]).collect()

    VARLOCIRAPTOR_PREPROCESS ( input.join( VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json), fasta, fai )

    scenario = file("./test.yml", checkIfExists: true)
    VARLOCIRAPTOR_CALLVARIANTS ( VARLOCIRAPTOR_PREPROCESS.out.vcf_gz, scenario )
}
