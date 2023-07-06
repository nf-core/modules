#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGNMENT } from '../../../../subworkflows/nf-core/alignment/main.nf'

workflow test_alignment {
    // load test data
    def bashScriptFile = new File('install_data.sh')

    def processBuilder = new ProcessBuilder('bash', bashScriptFile.toString())
    processBuilder.redirectOutput(ProcessBuilder.Redirect.INHERIT)
    processBuilder.redirectError(ProcessBuilder.Redirect.INHERIT)

    def process = processBuilder.start()
    process.waitFor()

    // channels enable parralle: https://www.nextflow.io/docs/latest/faq.html?highlight=parallel
    // test data 
    fastqs = [
    [[id:'gene', single_end:false], [params.test_data_msk['uncollapsed_bam_generation']['merged_fastq']['merged_1'], params.test_data_msk['uncollapsed_bam_generation']['merged_fastq']['merged_2']]]
    ]
    reference = [
        [id:'reference'], 
        file('test_nucleo/reference/')
    ]
    fastqs = ch_fastq = Channel.fromList(fastqs)

    // workflow 
    ALIGNMENT ( fastqs, reference, 1)
}
