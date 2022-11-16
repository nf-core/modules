#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// optional ags go below (if needed): https://nf-co.re/docs/contributing/tutorials/dsl2_modules_tutorial#passing-options-args
include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf'

//def TEST_IMG = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-0{6,7}.ome.tiff"
def TEST_IMG = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff"

workflow test_ashlar {

    println 'running workflow test_ashlar'

    /*
    // TO DO: test with input sample sheet like the following
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)=
    ASHLAR ( input )
    */

    input_channel = Channel.fromPath( TEST_IMG )
    //input_channel = Channel.fromPath( params.file_in )

    //ASHLAR ( file('/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff') )

    input =  [ [ id:'test' ],
               file(TEST_IMG, checkIfExists: true) ]

    ASHLAR ( input )
}
