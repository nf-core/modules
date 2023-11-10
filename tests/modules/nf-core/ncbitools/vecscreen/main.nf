#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB   } from '../../../../../modules/nf-core/blast/makeblastdb/main.nf'
include { NCBITOOLS_VECSCREEN } from '../../../../../modules/nf-core/ncbitools/vecscreen/main.nf'

workflow test_ncbitools_vecscreen {
    input = [ file(params.test_data['homo_sapiens']['pacbio']['primers'], checkIfExists: true) ]
    BLAST_MAKEBLASTDB   (  [[id:'test'], input] )
    NCBITOOLS_VECSCREEN (  [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db)
}
