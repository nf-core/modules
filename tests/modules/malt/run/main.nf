#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP } from '../../../../modules/unzip/main.nf' addParams( options: [:] )
include { MALT_BUILD } from  '../../../../modules/malt/build/main.nf' addParams( options: [:] )
include { MALT_RUN } from '../../../../modules/malt/run/main.nf' addParams( options: [:] )

workflow test_malt_run {

    fastas = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    seq_type = "DNA"
    map_db = file("https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jan2021.db.zip", checkIfExists: true)
    input = file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    mode = "BlastN"

    UNZIP ( map_db )
    MALT_BUILD ( fastas, seq_type, gff, UNZIP.out.result )
    MALT_RUN ( input, mode, MALT_BUILD.out.index )
}
