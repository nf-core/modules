#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SURVIVOR_SIMSV } from '../../../../../modules/nf-core/survivor/simsv/main.nf'

workflow test_survivor_simsv_parameters {
    
    SURVIVOR_SIMSV(
        [[],[]],
        [[],[]],
        [[],[]],
        [],
        []
    )
}

workflow test_survivor_simsv_vcf {

    fasta = [
        [id:"fasta"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ]

    fai = [
        [id:"fai"],
        file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true),
    ]

    parameters = Channel.of(
        "PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!\n" +
        "DUPLICATION_minimum_length: 100\n" +
        "DUPLICATION_maximum_length: 10000\n" +
        "DUPLICATION_number: 3\n" +
        "INDEL_minimum_length: 20\n" +
        "INDEL_maximum_length: 500\n" +
        "INDEL_number: 1\n" +
        "TRANSLOCATION_minimum_length: 1000\n" +
        "TRANSLOCATION_maximum_length: 3000\n" +
        "TRANSLOCATION_number: 0\n" +
        "INVERSION_minimum_length: 600\n" +
        "INVERSION_maximum_length: 800\n" +
        "INVERSION_number: 4\n" +
        "INV_del_minimum_length: 600\n" +
        "INV_del_maximum_length: 800\n" +
        "INV_del_number: 2\n" +
        "INV_dup_minimum_length: 600\n" +
        "INV_dup_maximum_length: 800\n" +
        "INV_dup_number: 2"
    ).collectFile(name:"parameters.txt", newLine:true)
    .map { [[id:"parameters"], it] }

    SURVIVOR_SIMSV(
        fasta,
        fai,
        parameters,
        1,
        1
    )
}

// A test with both outputs sadly isn't possible
// The default parameters contains translocation creation, which needs more than one chromosome
// The test date consists of only one chromosome


