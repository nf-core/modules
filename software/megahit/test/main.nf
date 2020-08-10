#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MEGAHIT } from '../main.nf'

dummy_file = file('DUMMY')


/*
 * Test with paired-end reads only
 */
workflow test_paired_end {
    reads = [[id: 'test', single_end: false, interleaved: false],
             ["${baseDir}/input/reads_R1.fastq.gz"],
             ["${baseDir}/input/reads_R2.fastq.gz"]]
                
    MEGAHIT (reads, [publish_dir:'test_paired_end'])
}

/*
 * Test with single-end reads only
 */
workflow test_single_end {

    reads = [[id: 'test', single_end: true, interleaved: false],
             ["${baseDir}/input/reads_single.fastq.gz"],
             [dummy_file]]
    MEGAHIT (reads, [ publish_dir:'test_single_end' ])
}

/*
 * Test with interleaved reads only
 */
workflow test_interleaved {

    reads = [[id: 'test', single_end: false, interleaved: true],
             ["${baseDir}/input/reads_interleaved.fastq.gz"],
             [dummy_file]]

    MEGAHIT (reads, [ publish_dir:'test_interleaved' ])
}

workflow {
    test_paired_end()
    test_single_end()
    test_interleaved()
}
