#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_RUN } from '../../../../../modules/nf-core/purecn/run/main.nf'

process STUB_PURECN_RUN {
    output:
    path("*interval_file.txt")    , emit: intervals
    path("*coverage.txt")         , emit: coverage
    path("*.vcf.gz")              , emit: vcf
    path("*normal_db.rds")        , emit: normal_db

    stub:
    """
    touch interval_file.txt
    touch coverage.txt
    touch test.vcf.gz
    touch normal_db.rds
    """
}

workflow test_purecn_run {
    
    input = [
        [ id:'test'],
        STUB_PURECN_RUN.out.intervals,
        STUB_PURECN_RUN.out.coverage,
        STUB_PURECN_RUN.out.vcf
    ]

    normal_db = STUB_PURECN_RUN.out.normal_db
    genome = "hg38"

    PURECN_RUN ( input, normal_db, genome )
}
