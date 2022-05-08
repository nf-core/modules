#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNTAR } from '../../../../modules/untar/main.nf'
include { SRATOOLS_FASTERQDUMP } from '../../../../modules/sratools/fasterqdump/main.nf'

workflow test_sratools_fasterqdump_single_end_with_input {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    archive = [ [], file(params.test_data['sarscov2']['illumina']['SRR13255544_tar_gz'], checkIfExists: true) ]
    UNTAR ( archive )

    def input = Channel.of([ id:'test_single_end', single_end:true ])
        .combine(UNTAR.out.untar.map{ it[1] })

    SRATOOLS_FASTERQDUMP(input, settings)
}

workflow test_sratools_fasterqdump_paired_end_with_input {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    archive = [ [], file(params.test_data['sarscov2']['illumina']['SRR11140744_tar_gz'], checkIfExists: true) ]
    UNTAR ( archive )

    def input = Channel.of([ id:'test_paired_end', single_end:false ])
        .combine(UNTAR.out.untar.map{ it[1] })

    SRATOOLS_FASTERQDUMP(input, settings)
}

workflow test_sratools_fasterqdump_single_end_without_input {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    archive = [ [], file(params.test_data['sarscov2']['illumina']['SRR13255544_tar_gz'], checkIfExists: true) ]
    UNTAR ( archive )

    def input = Channel.of([ id:'test_single_end', single_end:true ])
        .combine(UNTAR.out.untar.map{ it[1] })

    SRATOOLS_FASTERQDUMP(input, file('EXISTS'))
}

workflow test_sratools_fasterqdump_paired_end_without_input {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    archive = [ [], file(params.test_data['sarscov2']['illumina']['SRR11140744_tar_gz'], checkIfExists: true) ]
    UNTAR ( archive )

    def input = Channel.of([ id:'test_paired_end', single_end:false ])
        .combine(UNTAR.out.untar.map{ it[1] })

    SRATOOLS_FASTERQDUMP(input, file('EXISTS'))
}
