#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROTEUS_READPROTEINGROUPS } from '../../../../../modules/nf-core/proteus/readproteingroups/main.nf'

workflow test_proteus_readproteingroups {

    input = [
        [ id:'test' ], // meta map
        //file(params.test_data['proteomics']['maxquant']['samplesheet'], checkIfExists: true),    // samplesheet
        //file(params.test_data['proteomics']['maxquant']['proteingroups'], checkIfExists: true)  // intensities
        
    file('./delete_me/MaxQuant_samplesheet.tsv', checkIfExists: true),
    file('./delete_me/MaxQuant_proteinGroups.txt', checkIfExists: true)
    
    ]
    //ch_contrasts_file = Channel.fromPath(file(params.test_data['proteomics']['maxquant']['contrasts'], checkIfExists: true))    
    ch_contrasts_file = Channel.fromPath(file('./delete_me/MaxQuant_contrasts.csv', checkIfExists: true))
    ch_contrasts = ch_contrasts_file
        .map{it[1]}
        .splitCsv ( header:true, sep:',' )
        .map{
            tuple(
                [ id:'test' ],  // meta map
                it.variable     // contrast variable
            )
        }

    PROTEUS_READPROTEINGROUPS (
        input,
        ch_contrasts
    )
}
