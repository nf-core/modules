#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GFASTATS } from '../../../../modules/nf-core/gfastats/main.nf'

workflow test_gfastats_fasta_gfa {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    GFASTATS (
        input,
        'gfa',   // GFA output format
        '',      // No genome size
        '',      // No target
        [],      // No agp file
        [],      // No include bed file
        [],      // No exclude bed file
        []       // No swiss army knife instructions
    )
}

workflow test_gfastats_include_bed {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    bed   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    GFASTATS (
        input,
        'fastq', // Fastq output format
        '',      // No genome size
        '',      // No target
        [],      // No agp file
        bed,     // Include bed file
        [],      // No exclude bed file
        []       // No swiss army knife instructions
    )
}

workflow test_gfastats_exclude_bed {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    bed   = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

    GFASTATS (
        input,
        'fasta', // Fasta output format
        '',      // No genome size
        '',      // No target
        [],      // No agp file
        [],      // No include bed file
        bed,     // Exclude bed file
        []       // No swiss army knife instructions
    )
}

workflow test_gfastats_sak {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    sak   = Channel.of("""\
    RVCP\tMT192765.1
    """.stripIndent()).collectFile( name: 'assembly.sak' )

    GFASTATS (
        input,
        'fasta',      // Fasta output
        29800,        // Genome size
        'MT192765.1', // target
        [],           // No agp file
        [],           // No include bed file
        [],           // No exclude bed file
        sak           // Swiss army knife instructions
    )
}
