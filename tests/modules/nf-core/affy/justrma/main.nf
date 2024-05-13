#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AFFY_JUSTRMA } from '../../../../../modules/nf-core/affy/justrma/main.nf'
include { UNTAR        } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_affy_justrma {

    samples = file(params.test_data['homo_sapiens']['genome']['affy_array_samplesheet'], checkIfExists: true)
    cel_archive = file(params.test_data['homo_sapiens']['genome']['affy_array_celfiles_tar'], checkIfExists: true)
    
    def meta = [ id:'test' ]
    ch_samplesheet = Channel.of([ meta, samples ])
    ch_celfiles_archive = Channel.of([ meta, cel_archive ])

    UNTAR ( ch_celfiles_archive )

    ch_input = ch_samplesheet.join(UNTAR.out.untar)

    AFFY_JUSTRMA ( 
        ch_input,
        [[],[]] 
    )
}
