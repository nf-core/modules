#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHEWBBACA_ALLELECALL } from '../../../../../modules/nf-core/chewbbaca/allelecall/main.nf'
include { CHEWBBACA_CREATESCHEMA } from '../../../../../modules/nf-core/chewbbaca/createschema/main.nf'

workflow test_chewbbaca_allelecall {

    schema_input = [
        [ id:'test', single_end:false ], // meta map
        file("params.test_data['sarscov2']['illumina']['contigs_fasta']", checkIfExists: true)]
    ptf = []
    cds = []

    CHEWBBACA_CREATESCHEMA ( schema_input, ptf, cds )

    input = [
        [ id:'test', single_end:false ], // meta map
        file("/mnt/cidgoh-object-storage/hackathon/seqqc/genome_asm/GCA_010673125.1/GCA_010673125.1.fa", checkIfExists: true)
    ]
    schema = "/scratch/group_share/tmp/enterobase_senterica_cgmlst/"
    CHEWBBACA_ALLELECALL ( input, CHEWBBACA_CREATESCHEMA.out.schema )
}
