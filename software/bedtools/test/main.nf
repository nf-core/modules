#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

include { BEDTOOLS_SLOP } from '../slop/main.nf' addParams( options: [publish_dir:'test_bed_file'])
include { BEDTOOLS_REMOVEGENES } from '../removegenes/main.nf' addParams( options: [publish_dir:'test_bed_file'])
include { BEDTOOLS_ERNAS } from '../ernas/main.nf' addParams( options: [publish_dir:'test_bed_file'])
include { BEDTOOLS_TESTERNA } from '../testerna/main.nf' addParams( options: [publish_dir:'test_bed_file'])
include { BEDTOOLS_HISTONESTOBED} from '../histonestobed/main.nf' addParams( options: [publish_dir:'test_histones_to_bed'])
include { BEDTOOLS_ERNAGENEGROUPS} from '../ernagenegroups/main.nf' addParams( options: [publish_dir:'test_erna_gene_groups'])

// Define input channels

// Run the workflow
workflow test_bed_file {
    def input = []
    input = [ [ id:'test', single_end:true ], 
              [ file("${baseDir}/input/A.bed", checkIfExists: true),] ]

    BEDTOOLS_SLOP(
        input, 
        file("${baseDir}/input/genome.sizes", checkIfExists: true)
    )
    BEDTOOLS_REMOVEGENES(
        BEDTOOLS_SLOP.out.slopbed,
        file("${baseDir}/input/B.metatranscripts", checkIfExists: true)
        )
        
    BEDTOOLS_ERNAS(
        BEDTOOLS_REMOVEGENES.out.nogenesbed,
        file("${baseDir}/input/H3K27ac.bed", checkIfExists: true),
        file("${baseDir}/input/H3K4me1.bed", checkIfExists: true)
    )

    BEDTOOLS_TESTERNA(
        BEDTOOLS_ERNAS.out.ernabed, 
        file("${baseDir}/input/B.bed", checkIfExists: true)
    )

}

workflow test_histones_to_bed {
    def input = []
    input = [ [ id:'test', single_end:true ], 
              [ file("${baseDir}/input/test.single_end.sorted.bam", checkIfExists: true),] ]
    BEDTOOLS_HISTONESTOBED( input )
}


workflow {
    test_bed_file()
    test_histones_to_bed()
}


