#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Local default params
params.internal_process_name = 'umitools_dedup'

// Process definition
process umitools_dedup {
    publishDir "${params.outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-umitools:latest'

    input:
        tuple val(sample_id), path(bai), path(bam)
       
    output:
        tuple val(sample_id), path("*.dedup.bam"), emit: dedupBam
        path "*.dedup.log", emit: report

    script:

    // Init
    internal_prog = "umi_tools dedup"
    args = "--log=${sample_id}.dedup.log"

    // Check main args string exists and strip whitespace
    if(params.umitools_dedup_args) {
        ext_args = params.umitools_dedup_args
        internal_args += " " + ext_args.trim()
    }

    // Contruct CL line
    command = "${internal_prog} ${internal_args} -I $bam -S ${sample_id}.dedup.bam --output-stats=${sample_id}"

    // Log
    if (params.verbose){
        println ("[MODULE] umi_tools/dedup exec: " + internal_cl)
    }

    //SHELL
    """
    ${internal_cl}
    """
}
