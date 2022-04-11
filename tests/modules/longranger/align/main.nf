#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = "./"

include { LONGRANGER_ALIGN } from '../../../../modules/longranger/align/main.nf'

process DOWNLOAD_READS {
    script:
    """
    mkdir -p ${workDir}/10x
    for f in pEimTen1_S10_L007_R1_001.fastq.gz pEimTen1_S10_L007_R2_001.fastq.gz; do
        if [ ! -f  ${workDir}/10x/\$f ]
        then
            wget -P ${workDir}/10x https://darwin.cog.sanger.ac.uk/longranger_nf_test/10x/\$f
        fi
    done
    """
}

process DOWNLOAD_REF {
    script:
    """
    cd ${workDir}/
    wget https://darwin.cog.sanger.ac.uk/longranger_nf_test/refdata-pEimTen1.contigs.tar.gz
    tar -xzvf refdata-pEimTen1.contigs.tar.gz
    """
}

workflow test_longranger_align {
    DOWNLOAD_READS ()
    DOWNLOAD_REF ()

    LONGRANGER_ALIGN (
        [ [id : "pEimTen1"], file("${workDir}/10x/")],
        file("${workDir}/refdata-pEimTen1.contigs/"),
    )
}
