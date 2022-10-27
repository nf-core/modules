#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { GUNZIP        } from "$moduleDir/modules/nf-core/gunzip/main.nf"
include { GFFREAD       } from "$moduleDir/modules/nf-core/gffread/main.nf"
include { ULTRA_INDEX   } from "$moduleDir/modules/nf-core/ultra/index/main.nf"
include { ULTRA_ALIGN   } from "$moduleDir/modules/nf-core/ultra/align/main.nf"

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
