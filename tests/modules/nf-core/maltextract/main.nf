#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UNZIP as UNZIP_MALT        } from '../../../modules/unzip/main.nf'
include { UNZIP as UNZIP_MALTEXTRACT }  from '../../../modules/unzip/main.nf'
include { MALT_BUILD                 } from  '../../../modules/malt/build/main.nf'
include { MALT_RUN                   } from '../../../modules/malt/run/main.nf'
include { MALTEXTRACT                } from '../../../modules/maltextract/main.nf'

workflow test_maltextract {

    fastas = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    gff = []
    seq_type = "DNA"
    map_db = [ [], file("https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Jan2021.db.zip", checkIfExists: true) ]
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
    ]
    mode = "BlastN"
    taxon_list = file(params.test_data['sarscov2']['genome']['taxon_list_txt'], checkIfExists: true)
    ncbi_dir = [ [], file(params.test_data['sarscov2']['genome']['ncbi_taxmap_zip'], checkIfExists: true) ]

    UNZIP_MALT ( map_db )
    UNZIP_MALTEXTRACT ( ncbi_dir )
    MALT_BUILD ( fastas, seq_type, gff, UNZIP_MALT.out.unzipped_archive.map{ it[1] } )
    MALT_RUN ( input, mode, MALT_BUILD.out.index )
    ch_input_to_maltextract = MALT_RUN.out.rma6.map{ it[1] }
    MALTEXTRACT ( ch_input_to_maltextract, taxon_list, UNZIP_MALTEXTRACT.out.unzipped_archive.map{ it[1] })
}
