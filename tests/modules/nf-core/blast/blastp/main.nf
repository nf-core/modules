#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLAST_MAKEBLASTDB } from '../../../../../modules/nf-core/blast/makeblastdb/main.nf'
include { BLAST_BLASTP } from '../../../../../modules/nf-core/blast/blastp/main.nf'

workflow test_blast_blastp {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( [ [id:'test2'], input ] )
    out_ext = '' // empty test case to check default
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}

workflow test_blast_blastp_xml {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( [ [id:'test2'], input ] )
    out_ext = 'xml'
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}

workflow test_blast_blastp_tsv {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( [ [id:'test2'], input ] )
    out_ext = 'tsv'
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}

workflow test_blast_blastp_csv {

    input = [ file(params.test_data['sarscov2']['genome']['proteome_fasta'], checkIfExists: true) ]

    BLAST_MAKEBLASTDB ( [ [id:'test2'], input ] )
    out_ext = 'csv'
    BLAST_BLASTP ( [ [id:'test'], input ], BLAST_MAKEBLASTDB.out.db, out_ext )
}
