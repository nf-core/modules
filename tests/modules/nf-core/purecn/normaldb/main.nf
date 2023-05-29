#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_NORMALDB } from '../../../../../modules/nf-core/purecn/normaldb/main.nf'

workflow test_purecn_normaldb {

    input  = [
        [ id:'test' ],
        [file('https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/example_normal.txt.gz', checkIfExists: true),
        file('https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/example_normal2.txt.gz', checkIfExists: true)],
        [], []
    ]
    genome = 'hg38'
    assay  = 'illumina'

    PURECN_NORMALDB ( input, genome, assay )
}

workflow test_purecn_normaldb_normalvcf {

    input  = [
        [ id:'test' ],
        [file('https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/example_normal.txt.gz', checkIfExists: true),
        file('https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/example_normal2.txt.gz', checkIfExists: true)],
        [vcf: file('https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/normalpanel.vcf.gz', checkIfExists: true),
        tbi: file('https://raw.githubusercontent.com/lima1/PureCN/master/inst/extdata/normalpanel.vcf.gz.tbi', checkIfExists: true)]
    ]
    genome = 'hg38'
    assay  = 'illumina'

    PURECN_NORMALDB ( input, genome, assay )
}