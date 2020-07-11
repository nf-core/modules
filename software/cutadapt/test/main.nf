#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
include '../main.nf' params(params)

// Define input channels
input_fastqs = Channel.fromFilePairs('../../../test-datasets/tools/cutadapt/input/*_{1,2}.fastq' )

if(params.single_end){
    input_fastqs = Channel.from([
            ['SRR4238351', '../../../test-datasets/tools/cutadapt/input/SRR4238351_subsamp.fastq.gz'],
            ['SRR4238355', '../../../test-datasets/tools/cutadapt/input/SRR4238355_subsamp.fastq.gz'],
            ['SRR4238359', '../../../test-datasets/tools/cutadapt/input/SRR4238359_subsamp.fastq.gz'],
            ['SRR4238379', '../../../test-datasets/tools/cutadapt/input/SRR4238379_subsamp.fastq.gz']
        ]).map { row -> [ row[0], [ file(row[1]) ] ] }
}

// Run the workflow
workflow {
    cutadapt(input_fastqs)
}
