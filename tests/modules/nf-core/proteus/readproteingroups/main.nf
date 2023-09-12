#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROTEUS_READPROTEINGROUPS } from '../../../../../modules/nf-core/proteus/readproteingroups/main.nf'

workflow test_proteus_readproteingroups {

    input = [
        [ id:'test' ],                                                                              // meta map
        file(params.test_data['proteomics']['maxquant']['mq_samplesheet'], checkIfExists: true),    // samplesheet
        file(params.test_data['proteomics']['maxquant']['mq_proteingroups'], checkIfExists: true)   // intensities
    ]
    ch_contrasts_file = Channel.fromPath(file(params.test_data['proteomics']['maxquant']['mq_contrasts'], checkIfExists: true))    
    ch_contrasts = ch_contrasts_file
        .splitCsv(header:true, sep:',')
        .map{
            tuple(
                it,             // meta map
                it.variable     // contrast variable
            )
        }

    PROTEUS_READPROTEINGROUPS (
        input,
        ch_contrasts
    )

}
