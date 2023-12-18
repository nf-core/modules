#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { OATK } from '../../../../modules/nf-core/oatk/main.nf'

workflow test_oatk_mito {

    data_reads = Channel.of([[id:"ilDeiPorc1"],file(params.test_data['deilephila_porcellus']['mito']['hifi_reads'], checkIfExists: true)])
    hmm = file("https://raw.githubusercontent.com/c-zhou/OatkDB/main/v20230921/insecta_mito.fam", checkIfExists: true)
    hmm_h3f = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/insecta_mito.fam.h3f", checkIfExists: true)
    hmm_h3i = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/insecta_mito.fam.h3i", checkIfExists: true)
    hmm_h3m = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/insecta_mito.fam.h3m", checkIfExists: true)
    hmm_h3p = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/insecta_mito.fam.h3p", checkIfExists: true)

    OATK ( data_reads, hmm, hmm_h3f, hmm_h3i, hmm_h3m, hmm_h3p, [], [], [], [], [], [])
}

workflow test_oatk_pltd {

    data_reads = Channel.of([[id:"ddAraThal4"],file(params.test_data['arabidopsis_thaliana']['plastid']['hifi_reads'], checkIfExists: true)])  

    hmm = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/embryophyta_pltd.fam", checkIfExists: true)
    hmm_h3f = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/embryophyta_pltd.fam.h3f", checkIfExists: true)
    hmm_h3i = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/embryophyta_pltd.fam.h3i", checkIfExists: true)
    hmm_h3m = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/embryophyta_pltd.fam.h3m", checkIfExists: true)
    hmm_h3p = file("https://github.com/c-zhou/OatkDB/raw/main/v20230921/embryophyta_pltd.fam.h3p", checkIfExists: true)


    OATK ( data_reads, [], [], [], [], [], hmm, hmm_h3f, hmm_h3i, hmm_h3m, hmm_h3p, [])
}
