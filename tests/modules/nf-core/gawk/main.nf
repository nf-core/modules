#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAWK } from '../../../../modules/nf-core/gawk/main.nf'

workflow test_gawk {

    // A small example of a conversion of fasta index to bed format

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    GAWK (
        input,
        []
    )
}

workflow test_gawk_program_file {

    // A small example of a conversion of fasta index to bed format

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    program_file = Channel.of('BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}').collectFile(name:"program.txt")

    GAWK (
        input,
        program_file
    )
}
