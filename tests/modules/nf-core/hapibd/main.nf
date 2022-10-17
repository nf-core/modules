#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPIBD } from '../../../../modules/nf-core/hapibd/main.nf'

workflow test_hapibd {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hapibd/target.truth.vcf.gz", checkIfExists: true)
    ]
    map = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hapibd/target.map", checkIfExists: true)

    HAPIBD ( input, map, [] )
}

workflow test_hapibd_excludesamples {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hapibd/target.truth.vcf.gz", checkIfExists: true)
    ]
    map = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hapibd/target.map", checkIfExists: true)
    exclude = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hapibd/excludeSamples.txt", checkIfExists: true)

    HAPIBD ( input, map, exclude )
}
