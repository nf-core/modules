#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HLALA_TYPING } from   '../../../../../modules/nf-core/hlala/typing/main.nf'
include { SAMTOOLS_INDEX } from '../../../../../modules/nf-core/samtools/index/main.nf'

workflow test_hlala_typing {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]

    graph = Channel.fromPath(file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true))

    SAMTOOLS_INDEX(input)

    input_ch = Channel.from(input).toList()

    input_ch.concat(graph, SAMTOOLS_INDEX.out.bai)
                                            .flatten().toList()
                                            .map { meta, cram, graph, meta_1, bai ->
                                                    [meta, cram, bai, graph] }
                                            .set { hla_typing_input_ch }

    HLALA_TYPING ( hla_typing_input_ch )
}

 // extended tests running with actual test data (also see https://github.com/DiltheyLab/HLA-LA)
    
// Download files
// 
// wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
// wget https://www.dropbox.com/s/xr99u3vqaimk4vo/NA12878.mini.cram

// workflow test_hlala_typing {

//     input = [
//         [ id:'test' ], // meta map
//         file('/path/to/NA12878.mini.cram', checkIfExists: true)
//     ]

//     graph = Channel.fromPath("path/to/PRG_MHC_GRCh38_withIMGT")

//     SAMTOOLS_INDEX(input)

//     input_ch = Channel.from(input).toList()

//     input_ch.concat(graph, SAMTOOLS_INDEX.out.bai)
//                                             .flatten().toList()
//                                             .map { meta, cram, graph, meta_1, bai ->
//                                                     [meta, cram, bai, graph] }
//                                             .set { hla_typing_input_ch }

//     hla_typing_input_ch

//     HLALA_TYPING ( hla_typing_input_ch )
// }
