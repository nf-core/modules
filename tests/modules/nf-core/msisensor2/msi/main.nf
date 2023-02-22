#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSOR2_SCAN } from '../../../../../modules/nf-core/msisensor2/scan/main.nf'
include { MSISENSOR2_MSI } from '../../../../../modules/nf-core/msisensor2/msi/main.nf'

workflow test_msisensor2_msi_tumor_only {

    reference = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MSISENSOR2_SCAN ( reference, "outputfile" )

    input = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [],
        [],
        [],
    ]

    MSISENSOR2_MSI ( input, MSISENSOR2_SCAN.out.scan, [] )
}


workflow test_msisensor2_msi_tumor_normal {

    reference = [
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MSISENSOR2_SCAN ( reference, "outputfile" )

    input = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        [],
    ]

    MSISENSOR2_MSI ( input, MSISENSOR2_SCAN.out.scan, [] )
}

workflow test_msisensor2_msi_tumor_only_ml {

    input = [
        [ id:'test' ],
        file('https://github.com/niu-lab/msisensor2/raw/master/test/example.tumor.only.hg19.bam', checkIfExists: true),
        file('https://github.com/niu-lab/msisensor2/raw/master/test/example.tumor.only.hg19.bam.bai', checkIfExists: true),
        [],
        [],
        [],
    ]

    models = Channel.fromPath(
        [
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/016a16e12aca2bdba3713a3be76f72cd",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/02d42c2bda19aac304d6e86390c7f328",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/1030c0aa35ca5c263daeae866ad18632",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/15c3f5ec1c020d8f44283e40a2d9b6bb",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/15d6012f9a234b7adbbeecec524aea7d",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/2cf9a58f57e78b88acd86d792fe6a7b3",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/3ba04410c7ccbfc33e8b1b11d8132ae9",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/4431c9dc08be932c460a9e67192e7c57",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/4f5fa7bed97b48093375222d242fc982",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/71e6c0d59ea09d2a7acc566560841e34",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/8144b15900bba7086e86b31a0e1f8cfd",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/9bf6f7a544f369c3262a3a6f72cfdd7b",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/b8a36f2274b33cb0ed932e85cd1ddd5a",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/c08f164ded323a8c2606c408c555d73d",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/ceaa36ddbb76dc6eb6199ed946945788",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/e05d5da7208a924762311eddc4ec96c0",
            "https://github.com/niu-lab/msisensor2/raw/master/test/tmp/models_hg19_17sites/f8a20acf51ccb2b0ce6af42f24a8b5ef",
        ],
        checkIfExists: true
    )

    MSISENSOR2_MSI ( input, [], models.collect() )
}
