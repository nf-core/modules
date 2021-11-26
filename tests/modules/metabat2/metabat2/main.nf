#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { METABAT2_METABAT2 } from '../../../../modules/metabat2/metabat2/main.nf'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../../modules/metabat2/jgisummarizebamcontigdepths/main.nf'

workflow test_metabat2_no_depth {

    input_depth = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]

    Channel.fromPath(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        .map { it -> [[ id:'test', single_end:false ], it, []] }
        .set { input_metabat2 }
    
    METABAT2_METABAT2 ( input_metabat2 )
}

workflow test_metabat2_depth {

    input_depth = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( input_depth )

    Channel.fromPath(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
        .map { it -> [[ id:'test', single_end:false ], it] }
        .join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth)
        .set { input_metabat2 }

        METABAT2_METABAT2 ( input_metabat2 )
}
