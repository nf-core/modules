#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GALAH                       } from '../../../../modules/nf-core/galah/main.nf'
include { BIOAWK as BIOAWK_CHECKM     } from '../../../../modules/nf-core/bioawk/main.nf'
include { BIOAWK as BIOAWK_GENOMEINFO } from '../../../../modules/nf-core/bioawk/main.nf'
include { GUNZIP                      } from '../../../../modules/nf-core/gunzip/main.nf'


workflow test_galah {

    input = [
        [ id:'test' ], // meta map
        [file("https://github.com/nf-core/test-datasets/raw/magmap/testdata/GCA_002688505.1_ASM268850v1_genomic.fna.gz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/magmap/testdata/GCF_004296495.1_ASM429649v1_genomic.fna.gz", checkIfExists: true)]
    ]

    GALAH ( input, [], [] )

}
