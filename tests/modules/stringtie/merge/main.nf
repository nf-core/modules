#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STRINGTIE as STRINGTIE_FORWARD            } from '../../../../modules/stringtie/stringtie/main.nf'  addParams( options: [ publish_dir:'test_stringtie_forward' ] )
include { STRINGTIE as STRINGTIE_REVERSE            } from '../../../../modules/stringtie/stringtie/main.nf'  addParams( options: [ publish_dir:'test_stringtie_reverse' ] )
include { STRINGTIE_MERGE as STRINGTIE_FORWARD_MERGE} from '../../../../modules/stringtie/merge/main.nf'      addParams( options: [ publish_dir:'test_stringtie_forward_merge'] )
include { STRINGTIE_MERGE as STRINGTIE_REVERSE_MERGE} from '../../../../modules/stringtie/merge/main.nf'      addParams( options: [ publish_dir:'test_stringtie_reverse_merge'] )
/*
 * Test with forward strandedness
 */
workflow test_stringtie_forward_merge {
    input = [ [ id:'test', strandedness:'forward' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ] ]
    annotation_gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    
    STRINGTIE_FORWARD ( input, annotation_gtf )
    STRINGTIE_FORWARD.out.transcript_gtf
       .map { it -> it[1] }
       .set { stringtie_gtf }
    STRINGTIE_FORWARD_MERGE ( stringtie_gtf, annotation_gtf )
}

/*
 * Test with reverse strandedness
 */
workflow test_stringtie_reverse_merge {
    input = [ [ id:'test', strandedness:'reverse' ], // meta map
              [ file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true) ] 
            ]
    annotation_gtf   = file(params.test_data['homo_sapiens']['genome']['genome_gtf'], checkIfExists: true)
    
    STRINGTIE_REVERSE ( input, annotation_gtf )
    STRINGTIE_REVERSE.out.transcript_gtf
       .map { it -> it[1] }
       .set { stringtie_gtf }
    STRINGTIE_REVERSE_MERGE ( stringtie_gtf, annotation_gtf )
}
