#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf'

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
    ASHLAR ( input_channel )

    //ASHLAR ( file('/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff') )
}
