#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { COOLER_CLOAD }     from '../../../../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_BALANCE }   from '../../../../../modules/nf-core/cooler/balance/main.nf'
include { COOLTOOLS_PILEUP } from '../../../../../modules/nf-core/cooltools/pileup/main.nf'

workflow test_cooltools_pileup {

    input = [ [ id:'test', single_end:false ], // meta map
             file(params.test_data['generic']['cooler']['test_pairix_pair_gz'], checkIfExists: true),
             file(params.test_data['generic']['cooler']['test_pairix_pair_gz_px2'], checkIfExists: true),
             2000000
    ]

    sizes    = file(params.test_data['generic']['cooler']['hg19_chrom_sizes'], checkIfExists: true)

    COOLER_CLOAD ( input, sizes )
    COOLER_BALANCE ( COOLER_CLOAD.out.cool.map{[it[0], it[1], [:]]} )

    File bed  = new File("${workflow.workDir}/test.bed")
    bed.write("chr21\t0\t6000\r\nchr21\t10000\t320000\r\n")
    frag  = file(bed)

    COOLTOOLS_PILEUP ( COOLER_BALANCE.out.cool.combine([[:]]), frag )
}
