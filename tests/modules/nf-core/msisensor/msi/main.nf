#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MSISENSOR_SCAN } from '../../../../modules/msisensor/scan/main.nf'
include { MSISENSOR_MSI } from '../../../../modules/msisensor/msi/main.nf'

workflow test_msisensor_msi {

    def scaninput  = []
    scaninput = [ [ id:'test', single_end:false ], // meta map
                file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
    scan = MSISENSOR_SCAN ( scaninput )

    // IMPERFECT TEST:
    // USING SARS-COV2 DATA AS NORMAL:TUMOR PAIR THIS WILL SUFFICE TO
    // TEST MODULE EXECUTION, BUT NOT FUNCTIONALITY.
    // FUNCTIONALITY HAS BEEN TESTED MANUALY USING AUTHOR-PROVIDED TEST
    // DATA (https://github.com/ding-lab/msisensor/tree/master/test)
    def input      = []

    input = Channel.from([ [ id:'test', single_end:false ], // meta map
               file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam'],     checkIfExists: true),
               file(params.test_data['sarscov2']['illumina']['test_paired_end_methylated_sorted_bam_bai'], checkIfExists: true),
               file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'],         checkIfExists: true),
               file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'],     checkIfExists: true) ])
    //scan.txt.view()
    input.concat(scan.txt).collect().view()
    // BIT CLUMSY:
    MSISENSOR_MSI ( input.concat(scan.txt).collect() )
}
