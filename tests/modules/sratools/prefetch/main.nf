#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SRATOOLS_PREFETCH } from '../../../../modules/sratools/prefetch/main.nf'

workflow test_sratools_prefetch_with_settings_input {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    input = [
        [ id:'test', single_end:false ], // meta map
        'DRR000774'
    ]

    SRATOOLS_PREFETCH(input, settings)
}

workflow test_sratools_prefetch_without_settings_input {

    file(params.settings_path).mkdirs()
    def settings = file(params.settings_file)
    settings.text = "/LIBS/GUID = \"5b0d4b7d-88c7-4802-98fd-e3afd06feb32\"\n/libs/cloud/report_instance_identity = \"true\"\n"

    input = [
        [ id:'test', single_end:false ], // meta map
        'DRR000774'
    ]

    SRATOOLS_PREFETCH(input, file('EXISTS'))
}

