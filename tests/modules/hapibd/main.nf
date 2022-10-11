#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPIBD } from '../../../modules/hapibd/main.nf'

workflow test_hapibd {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/browning-lab/hap-ibd/master/test/target.truth.vcf.gz", checkIfExists: true)
    ]
    map = file("https://raw.githubusercontent.com/browning-lab/hap-ibd/master/test/target.map", checkIfExists: true)
    out = "hapIBDout"

    HAPIBD ( input, map, out, [] )
}

workflow test_hapibd_excludesamples {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://raw.githubusercontent.com/browning-lab/hap-ibd/master/test/target.truth.vcf.gz", checkIfExists: true)
    ]
    map = file("https://raw.githubusercontent.com/browning-lab/hap-ibd/master/test/target.map", checkIfExists: true)
    out = "hapIBDout"
    exclude = file("/Users/asma/Downloads/excludeSamples.txt", checkIfExists: true)

    HAPIBD ( input, map, out, exclude )
}
