#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { SOMATOSIM_SOMATOSIM } from '../../../../modules/nf-core/somatosim/somatosim'

workflow {
    // bam = Channel.fromPath('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.markduplicates.sorted.bam')
    
    bam = Channel.fromPath('/Users/stav/tools/SomatoSim/test_data/test_BAM.bam')
    bed = Channel.fromPath('/Users/stav/tools/SomatoSim/test_data/test_BED.bed')
    bam = bam.map{ it -> 
        meta = [:]
        meta.id = 'test_BAM'
        [meta, it]
    }

    
    SOMATOSIM_SOMATOSIM(bam, bed)
}