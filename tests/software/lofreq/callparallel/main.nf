#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX      } from '../../../../software/samtools/index/main.nf'      addParams( options: [:] )
include { SAMTOOLS_FAIDX      } from '../../../../software/samtools/faidx/main.nf'      addParams( options: [:] )
include { LOFREQ_CALLPARALLEL } from '../../../../software/lofreq/callparallel/main.nf' addParams( options: [:] )

workflow test_lofreq_callparallel {
    
    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ]
    
    fasta       = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    SAMTOOLS_INDEX ( input )
    SAMTOOLS_FAIDX ( fasta )

    LOFREQ_CALLPARALLEL ( input, SAMTOOLS_INDEX.out.bai, fasta, SAMTOOLS_FAIDX.out.fai )
}
