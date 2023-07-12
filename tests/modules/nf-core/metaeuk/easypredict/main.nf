#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MMSEQS_DATABASES                                  } from '../../../../../modules/nf-core/mmseqs/databases/main.nf'
include { METAEUK_EASYPREDICT as METAEUK_EASYPREDICT_FASTA  } from '../../../../../modules/nf-core/metaeuk/easypredict/main.nf'
include { METAEUK_EASYPREDICT as METAEUK_EASYPREDICT_MMSEQS } from '../../../../../modules/nf-core/metaeuk/easypredict/main.nf'

workflow test_metaeuk_easypredict_fasta {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fasta = [
        file(params.test_data['proteomics']['database']['yeast_ups'], checkIfExists: true)
    ]

    METAEUK_EASYPREDICT_FASTA ( input, fasta )

}

workflow test_metaeuk_easypredict_mmseqs {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    MMSEQS_DATABASES ( 'UniProtKB/Swiss-Prot' )
    database = MMSEQS_DATABASES.out.database

    METAEUK_EASYPREDICT_MMSEQS ( input, database )

}
