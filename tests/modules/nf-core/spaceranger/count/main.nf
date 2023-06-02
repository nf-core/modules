#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPACERANGER_COUNT } from '../../../../../modules/nf-core/spaceranger/count/main.nf'
include { UNTAR as DOWNLOAD_SPACERANGER_REFERENCE } from "../../../../../modules/nf-core/untar/main.nf"

workflow test_spaceranger_count {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    DOWNLOAD_SPACERANGER_REFERENCE(
        ['id': 'refdata-gex-GRCh38-2020-A'],
        file("https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz")
    )
    // SPACERANGER_COUNT ( input )
}
