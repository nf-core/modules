#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process umitools_dedup {
    publishDir "${params.outdir}/umitools/dedup",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-umitools:latest'

    input:
        tuple val(sample_id), path(bam)
       
    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedupBam
        tuple val(sample_id), path("${sample_id}.dedup.bam.bai"), emit: dedupBai
        path "*.dedup.log", emit: report

    script:

    // Init
    args = "--log=${sample_id}.dedup.log"

    // Check main args string exists and strip whitespace
    if(params.umitools_dedup_args) {
        ext_args = params.umitools_dedup_args
        args += " " + ext_args.trim()
    }

    // Contruct CL line
    dedup_command = "umi_tools dedup ${args} -I ${bam[0]} -S ${sample_id}.dedup.bam --output-stats=${sample_id}"

    // Log
    if (params.verbose){
        println ("[MODULE] umi_tools/dedup command: " + dedup_command)
    }

    //SHELL
    """
    ${dedup_command}
    samtools index ${sample_id}.dedup.bam
    """
}
