#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// import shovill
include {shovill} from '../main.nf' params(params)

// define input channel
readsPath = '../../../test-datasets/tools/shovill/input/SRR3609257_{1,2}.fastq.gz'
Channel
  .fromFilePairs( "${readsPath}", flat: true )
  .set{ ch_reads  }
 
// main workflow
workflow {
    shovill(ch_reads)
}
