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
        [['id': 'refdata-gex-GRCh38-2020-A'],
        file("https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz")]
    )
    ch_spaceranger_ref = DOWNLOAD_SPACERANGER_REFERENCE.out.untar.map({meta, ref -> ref})
    SPACERANGER_COUNT (
        [
            [
                id: "Visium_FFPE_Human_Ovarian_Cancer",
                slide: "V10L13-020",
                area: "D1"
            ],
            [
                file("./tests/modules/nf-core/spaceranger/count/testdata/human-ovarian-cancer-1-standard_v1_ffpe/Visium_FFPE_Human_Ovarian_Cancer_S1_L001_R1_001.fastq.gz"),
                file("./tests/modules/nf-core/spaceranger/count/testdata/human-ovarian-cancer-1-standard_v1_ffpe/Visium_FFPE_Human_Ovarian_Cancer_S1_L001_R2_001.fastq.gz")
            ],
            file("./tests/modules/nf-core/spaceranger/count/testdata/human-ovarian-cancer-1-standard_v1_ffpe/Visium_FFPE_Human_Ovarian_Cancer_image.jpg"),
            []
        ],
        ch_spaceranger_ref,
        []
    )
}
