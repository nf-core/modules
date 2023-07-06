#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_RUN } from '../../../../../modules/nf-core/purecn/run/main.nf'

process STUB_PURECN_RUN {
    output:
    path("*.txt")                , emit: intervals
    path("*.txt")                , emit: coverage
    path("*.vcf.gz")             , emit: vcf
    path("*.rds")                , emit: normal_db

    stub:
    """
    touch interval_file.txt
    touch coverage.txt
    touch test.vcf.gz
    touch normal_db.rds
    """
}

workflow test_purecn_run {

    STUB_PURECN_RUN()

    input = [
        [ id:'test'],
        file("interval_file.txt"),
        file("coverage.txt"),
        file("test.vcf.gz")
    ]

    normal_db = file("normal_db.rds")
    genome = "hg38"

    PURECN_RUN ( input, normal_db, genome )
}
