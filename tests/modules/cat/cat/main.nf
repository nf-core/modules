#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CAT_CAT } from '../../../../modules/cat/cat/main.nf' addParams( options: [:] )
include { CAT_CAT as CAT_CAT_SUFFIX } from '../../../../modules/cat/cat/main.nf' addParams( options: [suffix: ".fna"] )

workflow test_cat_ungzipped {
    
    input = [
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    ]

    CAT_CAT ( input )
}

workflow test_cat_gzipped {
    
    input = [
        file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true)
    ]

    CAT_CAT ( input )
}

workflow test_cat_ungzipped_fna {
    
    input = [
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    ]

    CAT_CAT_SUFFIX ( input )
}

workflow test_cat_gzipped_fna {
    
    input = [
        file(params.test_data['sarscov2']['genome']['genome_gff3_gz'], checkIfExists: true),
        file(params.test_data['sarscov2']['genome']['contigs_genome_maf_gz'], checkIfExists: true)
    ]

    CAT_CAT_SUFFIX ( input )
}
