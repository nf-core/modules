#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CAT_CAT } from '../../../../modules/cat/cat/main.nf'

workflow test_cat_unzipped_unzipped {

    input = [
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    ]

    CAT_CAT ( input, 'cat.txt' )
}

workflow test_cat_zipped_zipped {

    input = [
        file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true)
    ]

    CAT_CAT ( input, 'cat.txt.gz' )
}

workflow test_cat_zipped_unzipped {

    input = [
        file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true)
    ]

    CAT_CAT ( input, 'cat.txt' )
}

workflow test_cat_unzipped_zipped {

    input = [
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    ]

    CAT_CAT ( input, 'cat.txt.gz' )
}

workflow test_cat_one_file_unzipped_zipped {

    input = [
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
    ]

    CAT_CAT ( input, 'cat.txt.gz' )
}
