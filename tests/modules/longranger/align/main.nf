#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = "./"

include { LONGRANGER_ALIGN } from '../../../../modules/longranger/align/main.nf'

process DOWNLOAD_READS {
    script:
    """
    mkdir -p ${workDir}/10x
    wget -P ${workDir}/10x -c ${params.tol_test_data['test']['pEimTen1']['genomic_data']['tenx_r1_fastq_gz']}
    wget -P ${workDir}/10x -c ${params.tol_test_data['test']['pEimTen1']['genomic_data']['tenx_r2_fastq_gz']}
    """
}

process DOWNLOAD_REF {
    script:
    """
    cd ${workDir}/
    wget -c ${params.tol_test_data['small_genome']['pEimTen1']['assembly']['canu_contigs_longranger_mkref_targz']}
    tar -xzvf refdata-pEimTen1.contigs.tar.gz
    """
}

workflow test_longranger_align {
    DOWNLOAD_READS ()
    DOWNLOAD_REF ()

    LONGRANGER_ALIGN (
        [ [id : "pEimTen1"], file("${workDir}/refdata-pEimTen1.contigs/")],
        file("${workDir}/10x/"),
    )
}
