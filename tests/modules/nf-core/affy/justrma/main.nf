#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AFFY_JUSTRMA } from '../../../../../modules/nf-core/affy/justrma/main.nf'

workflow test_affy_justrma {

    def data_dir = 'https://raw.githubusercontent.com/pinin4fjords/test-datasets/affy/data/genomics/homo_sapiens/array_expression'

    samplesheet = file("${data_dir}/GSE38751.csv")
    sample_names = [ 'GSM949061_Dima_1', 'GSM949062_Dima_10', 'GSM949063_Dima_2', 'GSM949064_Dima_4', 'GSM949065_Dima_5', 'GSM949066_Dima_6', 'GSM949067_Dima_7', 'GSM949068_Dima_8', 'GSM949069_Dima_9', 'GSM949070_Dima_11' ]

    ch_input = Channel.from(sample_names)
        .map{file("${data_dir}/${it}.CEL.gz")}
        .toList()
        .map{
            [ [id:'test'], samplesheet, it ]
        }

    AFFY_JUSTRMA ( 
        ch_input,
        [[],[]] 
    )
}
