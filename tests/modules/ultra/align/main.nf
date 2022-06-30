#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP        } from '../../../../modules/gunzip/main.nf'
include { GFFREAD       } from '../../../../modules/gffread/main.nf'
include { ULTRA_INDEX   } from '../../../../modules/ultra/index/main.nf'
include { ULTRA_ALIGN   } from '../../../../modules/ultra/align/main.nf'

workflow test_ultra_align {

    input = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]

    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    gtf    = file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)

    GUNZIP ( input )
    GFFREAD ( gtf )
    ULTRA_INDEX ( genome, GFFREAD.out.gtf )
    ULTRA_ALIGN ( GUNZIP.out.gunzip, genome, ULTRA_INDEX.out.index )
}
