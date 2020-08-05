#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for umi_tools dedup...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.umitools_dedup_args = '--umi-separator=":"'
params.verbose = false

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include umitools_dedup from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testData = [
    ['sample1', "$baseDir/input/sample1.bam", "$baseDir/input/sample1.bai"],
    ['sample2', "$baseDir/input/sample2.bam", "$baseDir/input/sample2.bai"],
    ['sample3', "$baseDir/input/sample3.bam", "$baseDir/input/sample3.bai"],
    ['sample4', "$baseDir/input/sample4.bam", "$baseDir/input/sample4.bai"],
    ['sample5', "$baseDir/input/sample5.bam", "$baseDir/input/sample5.bai"],
    ['sample6', "$baseDir/input/sample6.bam", "$baseDir/input/sample6.bai"]
]

//Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_bam}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run dedup
    umitools_dedup ( ch_bam )
}

workflow.onComplete {
    def proc = "$baseDir/verify-checksum.sh $baseDir/../../../results/umitools/dedup/*.bam $baseDir/output/*.bam".execute()
    def b = new StringBuffer()
    proc.consumeProcessErrorStream(b)

    log.info proc.text

    errorString = b.toString()
    if(errorString != '')
        log.error errorString
        exit 1
}