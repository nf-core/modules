#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fastq_files = "out_L001_I1_001.fastq.gz out_L001_I2_001.fastq.gz out_L001_R1_001.fastq.gz out_L001_R2_001.fastq.gz"
params.read_structures = "8B 8B 150T 150T"

include { FQTK  } from '../../../../modules/nf-core/fqtk/main.nf'
include { UNTAR   } from '../../../../modules/nf-core/untar/main.nf'

workflow test_fqtk { 
    // Split params into lists and put into channels 
    ch_fastqs = Channel.from(params.fastq_files.split("\\s+"))
    ch_read_structures = Channel.from(params.read_structures.split("\\s+"))

    // Merge channel lists
    fastqs = ch_fastqs.merge( ch_read_structures ) 

    fastqs_with_paths = fastqs.combine(
        UNTAR ([
            [ id:'sim-data' ],
            file("https://github.com/nf-core/test-datasets/blob/demultiplex/testdata/sim-data/fastq.tar.gz?raw=true", checkIfExists: true)
        ]).untar.collect{it[1]}).toList()
    
    input = Channel.of ([ [ id:'sim-data'], // meta map
        file("https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/sim-data/fqtk_sample_metadata_subset.tsv", checkIfExists: true),
    ])

    ch_input = input.merge( fastqs_with_paths ) { a,b -> tuple(a[0], a[1], b)}

    FQTK ( ch_input )
}
