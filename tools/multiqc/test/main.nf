#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.outdir = "."
params.verbose = false
params.multiqc_args = ''

// include '../../../nf-core/module_testing/check_process_outputs.nf'
include '../main.nf'

if (params.verbose){
    println ("[WORKFLOW] MULTIQC ARGS: "          + params.multiqc_args)
}

multiqc_ch = Channel
  .fromPath( ['../../../test-datasets/*trimming_report.txt','../../../test-datasets/*fastqc.zip','../../../test-datasets/*screen.txt','../../../test-datasets/*bowtie2_stats.txt'] )
  .collect()   // collect() flattens all channels to single list
  // .view()   // view the files in the channel
  

// Run the workflow
workflow {

    main:
        // This is an example workflow for real reads aligned with Bowtie2. Just for illustration purposes
        
        // FASTQC          (file_ch, params.outdir, params.fastqc_args, params.verbose)
        // FASTQ_SCREEN    (file_ch, params.outdir, params.fastq_screen_args, params.verbose)
        // TRIM_GALORE     (file_ch, params.outdir, params.trim_galore_args, params.verbose)
        // FASTQC2         (TRIM_GALORE.out.reads, params.outdir, params.fastqc_args, params.verbose)
        // BOWTIE2         (TRIM_GALORE.out.reads, params.outdir, params.bowtie2_args, params.verbose)
     
        // merging channels for MultiQC
        // multiqc_ch = FASTQC.out.report.mix(
        //     TRIM_GALORE.out.report,
        //     FASTQ_SCREEN.out.report,
        //     FASTQC2.out.report,
        //     BOWTIE2.out.stats,
        // ).collect()     

        MULTIQC                          (multiqc_ch, params.outdir, params.multiqc_args, params.verbose)

        // .check_output() TODO
}