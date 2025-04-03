#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import the GCTA_GSMR module
include { GCTA_GSMR } from './'

// Define the workflow
workflow {
    // Create input channels
    exposure_ch = Channel.of([[id: "exposure"], [file(params.exposure)]])
    outcome_ch = Channel.of([[id: "outcome"], [file(params.outcome)]])


    // Run GCTA_GSMR

    GCTA_GSMR (
        exposure_ch,
        outcome_ch,
        params.reference
    )

    // Publish results
    GCTA_GSMR.out.gsmr_log.collectFile(name: 'gsmr_logs.txt', storeDir: params.outdir)
    GCTA_GSMR.out.gsmr_res.collectFile(name: 'gsmr_results.txt', storeDir: params.outdir)
    GCTA_GSMR.out.gsmr_effplot.collectFile(name: 'gsmr_effplots.txt', storeDir: params.outdir)
    GCTA_GSMR.out.versions.collectFile(name: 'software_versions.yml', storeDir: params.outdir)
}
