#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO as BUSCO_BACTE } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome --lineage_dataset bacteroidales_odb10'] )
include { BUSCO as BUSCO_CHR22 } from '../../../modules/busco/main.nf' addParams( options: [args: '--mode genome --offline'] )
include { UNTAR                } from '../../../modules/untar/main.nf' addParams( options: [:] )

// This tests genome decompression, empty input channels and data download
workflow test_busco_bacteroidales {
   input = [ [ id:'test' ], file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]    
   BUSCO_BACTE ( input,
                 [],
                 [] )
}

// This tests uncompressed genome, BUSCO lineage file provided via input channel, and offline mode
workflow test_busco_chr22 {
    input = [ [ id:'test' ],  file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true) ]
    lineage_dataset = [ file(params.test_data['homo_sapiens']['genome']['chr22_odb10_tar_gz'], checkIfExists: true) ]
    UNTAR(lineage_dataset)
    BUSCO_CHR22 ( input,
                  [],
                  UNTAR.out.untar )
}

