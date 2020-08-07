#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIMGALORE } from '../main.nf'

/*
 * Test with single-end data
 */
workflow test_single_end {

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${baseDir}/input/test_single_end.fastq.gz", checkIfExists: true) ] ]

    TRIMGALORE ( input, [ publish_dir:'test_single_end' ] )
}

/*
 * Test with paired-end data
 */
workflow test_paired_end {

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${baseDir}/input/test_R1.fastq.gz", checkIfExists: true),
                file("${baseDir}/input/test_R2.fastq.gz", checkIfExists: true) ] ]

    TRIMGALORE ( input, [ publish_dir:'test_paired_end' ] )
}

workflow {
    test_single_end()
    test_paired_end()
}




// #!/usr/bin/env nextflow
// nextflow.preview.dsl=2
//
// params.outdir = "."             // gets set in the nextflow.config files (to './results/trim_galore')
// params.verbose = false
// params.trim_galore_args = ''
// // trim_galore_args are best passed into the workflow in the following manner, e.g.:
// // --trim_galore_args="--clip_r1 10 --clip_r2 15 -j 2"
//
// if (params.verbose){
//     println ("[WORKFLOW] TRIM GALORE ARGS: "      + params.trim_galore_args)
// }
//
// // TODO: check the output files in some way
// // include '../../../tests/functions/check_process_outputs.nf'
// include '../main.nf'  // params (clip_r1: 6, clip_r2: 10) // how to pass additional parameters
//
// ch_read_files = Channel
//   .fromFilePairs('../../../test-datasets/test*{1,2}.fastq.gz',size:-1)
//   // .view()  // to check whether the input channel works
//
// workflow {
//
//     main:
//         TRIM_GALORE (ch_read_files, params.outdir, params.trim_galore_args, params.verbose)
//
// }
