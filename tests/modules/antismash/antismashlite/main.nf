#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GUNZIP } from '../../../modules/gunzip/main.nf'
include { UNTAR as UNTAR1 } from '../../../../modules/untar/main.nf'
include { UNTAR as UNTAR2 } from '../../../../modules/untar/main.nf'
include { UNTAR as UNTAR3 } from '../../../../modules/untar/main.nf'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../../modules/antismash/antismashlitedownloaddatabases/main.nf'
include { ANTISMASH_ANTISMASHLITE } from '../../../../modules/antismash/antismashlite/main.nf'
include { PROKKA } from '../../../modules/prokka/main.nf'

process STUB_ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES {
    input:
    path database_css
    path database_detection
    path database_modules

    output:
    path("antismash_db") , emit: database
    path("antismash_dir"), emit: antismash_dir

    stub:
    """
    mkdir antismash_db
    mkdir antismash_dir
    """
}

workflow test_antismashlite {
    input_genome = [
        [ id:'test' ], // meta map
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gbff_gz'], checkIfExists: true) ]

    input_antismash_db1 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz', checkIfExists: true)
    ]

    input_antismash_db2 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz', checkIfExists: true)
    ]

    input_antismash_db3 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz', checkIfExists: true)
    ]

    UNTAR1 ( input_antismash_db1 )
    UNTAR2 ( input_antismash_db2 )
    UNTAR3 ( input_antismash_db3 )
    ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] },  UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
    ANTISMASH_ANTISMASHLITE ( input_genome, ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.database, ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.antismash_dir )
}

workflow test_prokka_antismashlite {
    input_gunzip = [ [ id:'test' ], // meta map
            file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true) ]

    input_antismash_db1 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz', checkIfExists: true)
    ]

    input_antismash_db2 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz', checkIfExists: true)
    ]

    input_antismash_db3 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz', checkIfExists: true)
    ]

    GUNZIP ( input_gunzip )
    PROKKA ( GUNZIP.out.gunzip, [], [] )
    UNTAR1 ( input_antismash_db1 )
    UNTAR2 ( input_antismash_db2 )
    UNTAR3 ( input_antismash_db3 )
    STUB_ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] },  UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
    ANTISMASH_ANTISMASHLITE ( PROKKA.out.gbk, STUB_ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.database, STUB_ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.antismash_dir )
}
