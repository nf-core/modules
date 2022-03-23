#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.enable_conda = false

include { LONGRANGER_MKREF } from '../../../modules/longranger/mkref/main.nf' 
include { LONGRANGER_ALIGN } from '../../..//modules/longranger/align/main.nf' 

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

workflow {
    DOWNLOAD_READS()
    LONGRANGER_MKREF([ [id : "longranger test"], 
        file('https://darwin.cog.sanger.ac.uk/longranger_nf_test/pEimTen1.contigs.fasta', checkIfExists: true) ])
    reference = LONGRANGER_MKREF.out.folder
    LONGRANGER_ALIGN(
        [ [id : "longranger test"], "pEimTen1"], 
        [ [id : "longranger test"], file("${workDir}/10x/")], 
        reference)
}
