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
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.Aligned.sortedByCoord.out.bam"],
  ['Sample2', "$baseDir/input/prpf8_ctrl_rep2.Aligned.sortedByCoord.out.bam"],
  ['Sample3', "$baseDir/input/prpf8_ctrl_rep4.Aligned.sortedByCoord.out.bam"],
  ['Sample4', "$baseDir/input/prpf8_eif4a3_rep1.Aligned.sortedByCoord.out.bam"],
  ['Sample5', "$baseDir/input/prpf8_eif4a3_rep2.Aligned.sortedByCoord.out.bam"],
  ['Sample6', "$baseDir/input/prpf8_eif4a3_rep4.Aligned.sortedByCoord.out.bam"]
]

testMetaDataBai = [
  ['Sample1', "$baseDir/input/prpf8_ctrl_rep1.Aligned.sortedByCoord.out.bai"],
  ['Sample2', "$baseDir/input/prpf8_ctrl_rep2.Aligned.sortedByCoord.out.bai"],
  ['Sample3', "$baseDir/input/prpf8_ctrl_rep4.Aligned.sortedByCoord.out.bai"],
  ['Sample4', "$baseDir/input/prpf8_eif4a3_rep1.Aligned.sortedByCoord.out.bai"],
  ['Sample5', "$baseDir/input/prpf8_eif4a3_rep2.Aligned.sortedByCoord.out.bai"],
  ['Sample6', "$baseDir/input/prpf8_eif4a3_rep4.Aligned.sortedByCoord.out.bai"]
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