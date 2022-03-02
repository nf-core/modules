#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { METABAT2_METABAT2 } from '../../../../modules/metabat2/metabat2/main.nf'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../../modules/metabat2/jgisummarizebamcontigdepths/main.nf'
include { DASTOOL_SCAFFOLDS2BIN } from '../../../../modules/dastool/scaffolds2bin/main.nf'
include { DASTOOL_DASTOOL } from '../../../../modules/dastool/dastool/main.nf'

workflow test_dastool_dastool {

    input_depth = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['bacteroides_fragilis']['illumina']['test1_paired_end_sorted_bam'], checkIfExists: true),
              file(params.test_data['bacteroides_fragilis']['illumina']['test1_paired_end_sorted_bam_bai'], checkIfExists: true) ]

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( input_depth )

    Channel.fromPath(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)
        .map { it -> [[ id:'test', single_end:false ], it] }
        .join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth)
        .set { input_metabat2 }

    METABAT2_METABAT2 ( input_metabat2 )

    DASTOOL_SCAFFOLDS2BIN ( METABAT2_METABAT2.out.fasta.collect(), "fa")

    Channel.of([ [ id:'test', single_end:false ], // meta map
                      file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)])
        .join(DASTOOL_SCAFFOLDS2BIN.out.scaffolds2bin)
        .set {input_dastool}
    

    DASTOOL_DASTOOL ( input_dastool, [], [], [] )
}
