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
    DASTOOL_SCAFFOLDS2BIN ( MAXBIN2.out.binned_fastas.collect(), "binning_software_name", "fasta")


    def input_dastool = Channel.of([ [ id:'test', single_end:false ], // meta map
                      file(params.test_data['bacteroides_fragilis']['illumina']['test1_contigs_fa_gz'], checkIfExists: true)])
    def proteins_empty = Channel.of([ [ id:'test', single_end:false ], [] ])

    input_dastool
        .join(DASTOOL_SCAFFOLDS2BIN.out.scaffolds2bin)
        .join(proteins_empty)
        .view()
        .set {input_dastool_full}
    

    DASTOOL_DASTOOL ( input_dastool_full, "binning_software_name", [] )
}
