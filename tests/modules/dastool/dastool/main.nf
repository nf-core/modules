#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MAXBIN2 } from '../../../modules/maxbin2/main.nf' addParams( options: [:] )
include { DASTOOL_SCAFFOLDS2BIN } from '../../../../modules/dastool/scaffolds2bin/main.nf' addParams( options: [:] )
include { DASTOOL_DASTOOL } from '../../../../modules/dastool/dastool/main.nf' addParams( options: [args: '--score_threshold 0'] )

workflow test_dastool_dastool {

    input_maxbin = [ [ id:'test', single_end:false ], // meta map
                     file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true),
                     file(params.test_data['bacteroides_fragilis']['illumina']['test1_1_fastq_gz'], checkIfExists: true),
                     []
                   ]
    MAXBIN2 ( input_maxbin )
    DASTOOL_SCAFFOLDS2BIN ( MAXBIN2.out.binned_fastas.collect(), "fasta")

    Channel.of([ [ id:'test', single_end:false ], // meta map
                      file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)])
        .join(DASTOOL_SCAFFOLDS2BIN.out.scaffolds2bin)
        .set {input_dastool}
    

    DASTOOL_DASTOOL ( input_dastool, [], [], [] )
}
