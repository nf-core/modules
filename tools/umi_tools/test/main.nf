#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test for umi_tools dedup...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.umitools_dedup_args = '--umi-separator=":"'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include umitools_dedup from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Meta data
testMetaDataBam = [
  ['sample1', "$baseDir/input/sample1.bam"],
  ['sample2', "$baseDir/input/sample2.bam"],
  ['sample3', "$baseDir/input/sample3.bam"],
  ['sample4', "$baseDir/input/sample4.bam"],
  ['sample5', "$baseDir/input/sample5.bam"],
  ['sample6', "$baseDir/input/sample6.bam"]
]

testMetaDataBai = [
  ['sample1', "$baseDir/input/sample1.bai"],
  ['sample2', "$baseDir/input/sample2.bai"],
  ['sample3', "$baseDir/input/sample3.bai"],
  ['sample4', "$baseDir/input/sample4.bai"],
  ['sample5', "$baseDir/input/sample5.bai"],
  ['sample6', "$baseDir/input/sample6.bai"]
]

//Bam input channel
Channel
  .from(testMetaDataBam)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ]}
  .set {ch_test_meta_bam}

//BamBai input channel
Channel
  .from(testMetaDataBai)
  .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
  .join(ch_test_meta_bam)
  .set {ch_test_meta_bambai} 

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run dedup
    umitools_dedup ( ch_test_meta_bambai )

    // Collect file names and view output
    umitools_dedup.out.dedupBam | view
}