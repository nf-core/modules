#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ANTISMASH_ANTISMASHLITE } from '../../../../modules/antismash/antismashlite/main.nf'
include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../../modules/antismash/antismashlitedownloaddatabases/main.nf'
include { GUNZIP as GUNZIP1 } from '../../../../modules/gunzip/main.nf'
include { GUNZIP as GUNZIP2 } from '../../../../modules/gunzip/main.nf'
include { UNTAR as UNTAR1 } from '../../../../modules/untar/main.nf'
include { UNTAR as UNTAR2 } from '../../../../modules/untar/main.nf'
include { UNTAR as UNTAR3 } from '../../../../modules/untar/main.nf'

workflow test_antismashlite {
    genome_fna = [
        [ id:'test' ],
        file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
    ]

    genome_gff = [
        [],
        file(params.test_data['bacteroides_fragilis']['genome']['genome_gff_gz'], checkIfExists: true)
    ]

    antismash_db1 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz', checkIfExists: true)
    ]

    antismash_db2 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz', checkIfExists: true)
    ]

    antismash_db3 = [
        [],
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz', checkIfExists: true)
    ]

    // input_db = [
    //     file('/home/jasmin/antismash_db/', checkIfExists: true)
    // ]

    // input_dir = [
    //     file('/home/jasmin/antismash_dir/', checkIfExists: true)
    // ]

    GUNZIP1 ( genome_fna )
    GUNZIP2 ( genome_gff )
    UNTAR1 ( antismash_db1 )
    UNTAR2 ( antismash_db2 )
    UNTAR3 ( antismash_db3 )
    ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] },  UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
    ANTISMASH_ANTISMASHLITE ( GUNZIP1.out.gunzip, ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.database, ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES.out.antismash_dir, GUNZIP2.out.gunzip.map{ it[1] } )
    // ANTISMASH_ANTISMASHLITE ( GUNZIP1.out.gunzip, input_db, input_dir, GUNZIP2.out.gunzip.map{ it[1] } )

}
// #!/usr/bin/env nextflow

// nextflow.enable.dsl = 2

// // include { UNTAR as UNTAR1 } from '../../../../modules/untar/main.nf'
// // include { UNTAR as UNTAR2 } from '../../../../modules/untar/main.nf'
// // include { UNTAR as UNTAR3 } from '../../../../modules/untar/main.nf'
// // // include { ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES } from '../../../modules/antismash/antismashlitedownloaddatabases/main.nf'
// include { ANTISMASH_ANTISMASHLITE } from '../../../../modules/antismash/antismashlite/main.nf'

// workflow test_antismashlite {
//     input_genome = [
//         [ id:'test' ], // meta map
//         file('/home/jasmin/Downloads/BAN001-megahit.gbk', checkIfExists: true) ]

//     // input_antismash_db1 = [
//     //     [],
//     //     file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/css.tar.gz', checkIfExists: true)
//     // ]

//     // input_antismash_db2 = [
//     //     [],
//     //     file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/detection.tar.gz', checkIfExists: true)
//     // ]

//     // input_antismash_db3 = [
//     //     [],
//     //     file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/antismash/modules.tar.gz', checkIfExists: true)
//     // ]

//     input_db = [
//         [],
//         file('/home/jasmin/antismash_db/', checkIfExists: true)
//     ]

//     input_dir = [
//         [],
//         file('/home/jasmin/antismash_dir', checkIfExists: true)
//     ]

//     // input_gff = [
//     //     [],
//     //     file('/home/jasmin/Downloads/BAN001.gff', checkIfExists: true)
//     // ]

//     // UNTAR1 ( input_antismash_db1 )
//     // UNTAR2 ( input_antismash_db2 )
//     // UNTAR3 ( input_antismash_db3 )
//     // ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES ( UNTAR1.out.untar.map{ it[1] },  UNTAR2.out.untar.map{ it[1] }, UNTAR3.out.untar.map{ it[1] } )
//     ANTISMASH_ANTISMASHLITE ( input_genome, input_db, input_dir )
// }
