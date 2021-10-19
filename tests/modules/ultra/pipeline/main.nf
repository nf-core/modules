#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRA_PIPELINE } from '../../../../modules/ultra/pipeline/main.nf' addParams( options: [:] )

workflow test_ultra_pipeline {
    def wget1    = "wget -O " + workflow.workDir + "/genome.gtf "         + params.test_data['homo_sapiens']['genome']['genome_gtf']
    def wget2    = "wget -O " + workflow.workDir + "/test_hifi.fastq.gz " + params.test_data['homo_sapiens']['pacbio']['hifi']
    def sort_gtf = "sort -n -k 1,1 -k4,4 -k5,5 -o " + workflow.workDir + "/genome_sort.gtf " + workflow.workDir + "/genome.gtf"
    def gunzip   = "gunzip " + workflow.workDir + "/test_hifi.fastq.gz"

    wget1.execute().waitFor()
    wget2.execute().waitFor()
    sort_gtf.execute().waitFor()
    gunzip.execute().waitFor()

    input =
        [
            [ id:'test', single_end:false ], // meta map
            file(workflow.workDir + "/test_hifi.fastq", checkIfExists: true)
        ]

    gtf    = file(workflow.workDir + "/genome_sort.gtf",                      checkIfExists: true)
    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    ULTRA_PIPELINE ( input, genome, gtf )
}
