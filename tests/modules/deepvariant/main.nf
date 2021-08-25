#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEEPVARIANT } from '../../../modules/deepvariant/main.nf' addParams( options: ["args": "--regions=\"chr20:10,000,000-10,010,000\" --model_type=WGS"] )

workflow test_deepvariant {
    
    // input = [ [[ id:'test', single_end:false ], // meta map
    //            "https://storage.googleapis.com/deepvariant/quickstart-testdata/NA12878_S1.chr20.10_10p1mb.bam",
    //            "https://storage.googleapis.com/deepvariant/quickstart-testdata/NA12878_S1.chr20.10_10p1mb.bam.bai"],
    //          "https://storage.googleapis.com/deepvariant/quickstart-testdata/ucsc.hg19.chr20.unittest.fasta",
    //          "https://storage.googleapis.com/deepvariant/quickstart-testdata/ucsc.hg19.chr20.unittest.fasta.fai"
    // ]

    bam_tuple_ch = Channel.of([[ id:'test', single_end:false ], // meta map
               "https://storage.googleapis.com/deepvariant/quickstart-testdata/NA12878_S1.chr20.10_10p1mb.bam",
                               "https://storage.googleapis.com/deepvariant/quickstart-testdata/NA12878_S1.chr20.10_10p1mb.bam.bai"])

    fasta_tuple_ch = Channel.of(["https://storage.googleapis.com/deepvariant/quickstart-testdata/ucsc.hg19.chr20.unittest.fasta",
                  "https://storage.googleapis.com/deepvariant/quickstart-testdata/ucsc.hg19.chr20.unittest.fasta.fai"])

    DEEPVARIANT ( bam_tuple_ch, fasta_tuple_ch )
}
