#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CAT_CAT                        } from '../../../../modules/cat/cat/main.nf'
include { CAT_CAT as CAT_UNZIPPED_ZIPPED } from '../../../../modules/cat/cat/main.nf'
include { CAT_CAT as CAT_ZIPPED_UNZIPPED } from '../../../../modules/cat/cat/main.nf'

workflow test_cat_unzipped_unzipped {

    input = [
        [ id:'test', single_end:true ], // meta map
        [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true) ]
    ]

    CAT_CAT ( input )
}

workflow test_cat_zipped_zipped {

    input = [
         [ id:'test', single_end:true ], // meta map
        [file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true)]
    ]

    CAT_CAT ( input )
}

workflow test_cat_zipped_unzipped {

    input = [
        [ id:'test', single_end:true ], // meta map
        [file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true)]
    ]

    CAT_ZIPPED_UNZIPPED ( input )
}

workflow test_cat_unzipped_zipped {

    input = [
        [ id:'test', single_end:true ], // meta map
        [file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)]
    ]

    CAT_UNZIPPED_ZIPPED ( input )
}

workflow test_cat_one_file_unzipped_zipped {

    input = [
        [ id:'test', single_end:true ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    CAT_UNZIPPED_ZIPPED ( input )
}
