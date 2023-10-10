#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROTEUS_READPROTEINGROUPS } from '../../../../../modules/nf-core/proteus/readproteingroups/main.nf'

workflow test_proteus_readproteingroups {

    ch_input = Channel.of(
        [
            file(params.test_data['proteomics']['maxquant']['mq_samplesheet'], checkIfExists: true),    // samplesheet
            file(params.test_data['proteomics']['maxquant']['mq_proteingroups'], checkIfExists: true)   // intensities
        ]
    )
    
    ch_contrasts_file = Channel.fromPath(file(params.test_data['proteomics']['maxquant']['mq_contrasts'], checkIfExists: true))    
    ch_contrasts = ch_contrasts_file
        .splitCsv(header:true, sep:',')
        .map{
            tuple(
                [ id: it.variable ]                                                                     // metamap with contrast variable,
            )
        }
    ch_proteus_in = ch_contrasts.combine(ch_input)

    PROTEUS_READPROTEINGROUPS (
        ch_proteus_in
    )

}
