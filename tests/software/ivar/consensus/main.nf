#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.publish_dir_mode = 'link'
params.enable_conda = false

include { IVAR_CONSENSUS } from '../../../../software/ivar/consensus/main.nf' addParams([:])

workflow test_ivar_consensus {
    
    bam = file("${launchDir}/tests/data/bam/test-sc2-artic-v3.bam", checkIfExists: true)

    def input = []
    input = [ [ id:'test'], 
                file("${launchDir}/tests/data/bam/test-sc2-artic-v3-sorted-trimmed.bam", checkIfExists: true) ]
  IVAR_CONSENSUS ( input )
}
