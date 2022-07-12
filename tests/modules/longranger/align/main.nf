#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LONGRANGER_ALIGN } from '../../../../modules/longranger/align/main.nf'

process DOWNLOAD_READS {
    input:
    output:
    path("10x/"), emit: folder
    path("10x/*fastq.gz"), emit: reads
    script:
    """
    mkdir -p 10x
    wget -P 10x -c ${params.tol_test_data['test']['pEimTen1']['genomic_data']['tenx_r1_fastq_gz']}
    wget -P 10x -c ${params.tol_test_data['test']['pEimTen1']['genomic_data']['tenx_r2_fastq_gz']}
    """
}

process DOWNLOAD_REF {
    input:
    output:
    path("refdata-pEimTen1.contigs"), emit: contigs
    script:
    """
    wget -c ${params.tol_test_data['small_genome']['pEimTen1']['assembly']['canu_contigs_longranger_mkref_targz']}
    tar -xzvf refdata-pEimTen1.contigs.tar.gz
    """
}

workflow test_longranger_align {
    DOWNLOAD_READS ()
    DOWNLOAD_REF ()
    DOWNLOAD_READS.out.folder.view()
    DOWNLOAD_REF.out.contigs.view()
    LONGRANGER_ALIGN (
        DOWNLOAD_REF.out.contigs.map{ it -> [ [id : "pEimTen1"], it ]},
        DOWNLOAD_READS.out.folder,
    )
}
