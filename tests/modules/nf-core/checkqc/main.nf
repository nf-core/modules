#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKQC } from '../../../../modules/nf-core/checkqc/main.nf'
include { UNTAR } from '../../../../modules/nf-core/untar/main.nf'

workflow test_checkqc {

    run_dir_tar = [
        [],
        file('/home/matilda/workspace/test-datasets/data/genomics/homo_sapiens/illumina/bcl/flowcell_checkqc.tar.gz', checkIfExists: true)
    ]

    UNTAR ( run_dir_tar )

    CHECKQC ( UNTAR.out.untar.map{ it[1] }, [])
}
