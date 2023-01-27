#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AFFY_JUSTRMA } from '../../../../../modules/nf-core/affy/justrma/main.nf'
include { UNTAR        } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_affy_justrma {

    def data_dir = 'https://raw.githubusercontent.com/pinin4fjords/test-datasets/affy/data/genomics/homo_sapiens/array_expression'
    def meta = [ id:'test' ]
    
    samples = file("${data_dir}/GSE38751.csv")
    cel_archive = file("${data_dir}/GSE38751_RAW.tar")

    ch_samplesheet = Channel.of([ meta, samples ])
    ch_celfiles_archive = Channel.of([ meta, cel_archive ])

    UNTAR ( ch_celfiles_archive )

    ch_input = ch_samplesheet.join(UNTAR.out.untar)

    AFFY_JUSTRMA ( 
        ch_input,
        [[],[]] 
    )
}
