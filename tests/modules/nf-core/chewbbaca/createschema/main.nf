#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHEWBBACA_CREATESCHEMA } from '../../../../../modules/nf-core/chewbbaca/createschema/main.nf'

workflow test_chewbbaca_createschema {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true)]
    ptf = []
    cds = []

    CHEWBBACA_CREATESCHEMA ( input, ptf, cds)
}

workflow test_chewbbaca_createschema_multi {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true),
            file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)] ]
    ptf = []
    cds = []

    CHEWBBACA_CREATESCHEMA ( input, ptf, cds)
}

workflow test_chewbbaca_createschema_gz {

    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['sarscov2']['illumina']['contigs_fasta'], checkIfExists: true),
            file(params.test_data['sarscov2']['genome']['genome_fasta_gz'], checkIfExists: true)] ]
    ptf = []
    cds = []

    CHEWBBACA_CREATESCHEMA ( input, ptf, cds)
}
