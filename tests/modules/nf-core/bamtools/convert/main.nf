#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_EXT_ERROR }   from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_NOEXT_ERROR } from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_BED }         from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_FASTA }       from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_FASTQ }       from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_JSON }        from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_PILEUP }      from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_SAM }         from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"
include { BAMTOOLS_CONVERT as BAMTOOLS_CONVERT_YAML }        from "$moduleDir/modules/nf-core/bamtools/convert/main.nf"

workflow test_bamtools_convert_ext_error {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_EXT_ERROR ( input )
}

workflow test_bamtools_convert_noext_error {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_NOEXT_ERROR ( input )
}

workflow test_bamtools_convert_bed {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_BED ( input )
}

workflow test_bamtools_convert_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_FASTA ( input )
}

workflow test_bamtools_convert_fastq {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_FASTQ ( input )
}

workflow test_bamtools_convert_json {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_JSON ( input )
}

workflow test_bamtools_convert_pileup {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_PILEUP ( input )
}

workflow test_bamtools_convert_sam {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_SAM ( input )
}

workflow test_bamtools_convert_yaml {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    BAMTOOLS_CONVERT_YAML ( input )
}

