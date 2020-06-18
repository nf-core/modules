#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test for BAM deduplication")

/* Module inclusions 
--------------------------------------------------------------------------------------*/

include dedup from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testMetaDataBam = [
  ['Sample1', "$baseDir/input/sample1.bam"],
  ['Sample2', "$baseDir/input/sample2.bam"],
  ['Sample3', "$baseDir/input/sample3.bam"],
  ['Sample4', "$baseDir/input/sample4.bam"],
  ['Sample5', "$baseDir/input/sample5.bam"],
  ['Sample6', "$baseDir/input/sample6.bam"]
]

testMetaDataBai = [
  ['Sample1', "$baseDir/input/sample1.bai"],
  ['Sample2', "$baseDir/input/sample2.bai"],
  ['Sample3', "$baseDir/input/sample3.bai"],
  ['Sample4', "$baseDir/input/sample4.bai"],
  ['Sample5', "$baseDir/input/sample5.bai"],
  ['Sample6', "$baseDir/input/sample6.bai"]
]

// Create channels of test data 

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

// Run workflow
workflow {

    // Run dedup
    dedup ( ch_test_meta_bambai )

    // Collect file names and view output
    dedup.out.dedupBam | view
}