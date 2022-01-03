#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC                      } from '../../../../modules/fastqc/main.nf'
include { MULTIQC                     } from '../../../../modules/multiqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../../../modules/custom/dumpsoftwareversions/main.nf'

workflow fastqc1 {
    take:
    input

    main:
    FASTQC ( input )

    emit:
    versions = FASTQC.out.versions
}

workflow fastqc2 {
    take:
    input

    main:
    FASTQC ( input )

    emit:
    versions = FASTQC.out.versions
    zip = FASTQC.out.zip
}

workflow test_custom_dumpsoftwareversions {
    input = [
        [ id: 'test', single_end: false ],
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    // Using subworkflows to ensure that the script can properly handle
    // cases where subworkflows have a module with the same name.
    fastqc1 ( input )
    fastqc2 ( input )
    MULTIQC ( fastqc2.out.zip.collect { it[1] } )

    fastqc1
        .out
        .versions
        .mix(fastqc2.out.versions)
        .mix(MULTIQC.out.versions)
        .set { ch_software_versions }

    CUSTOM_DUMPSOFTWAREVERSIONS ( ch_software_versions.collectFile() )
}
