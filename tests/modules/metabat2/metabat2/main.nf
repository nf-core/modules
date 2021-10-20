#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METABAT2_METABAT2 } from '../../../../modules/metabat2/metabat2/main.nf' addParams( options: [args: '--minContig 1500 --minCV 0.1 --minCVSum 0.1 --minClsSize 10 --minS 2'] )
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../../modules/metabat2/jgisummarizebamcontigdepths/main.nf' addParams( options: [:] )

workflow test_metabat2_metabat2 {

    input_depth = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]

    

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( input_depth )

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .join(
            Channel.fromPath(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
                .map { it -> [[ id:'test', single_end:false ], it]}
        ). set { input_binning }

    METABAT2_METABAT2 ( input_binning )
}