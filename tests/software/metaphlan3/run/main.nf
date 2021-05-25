#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
//  test data locations in nf config (will tidy before final submission :) )
// To get around metaphlan auto-download, I used the following link (https://zenodo.org/record/4629921/files/metaphlan_databases.tar.gz) to download the tar file (helpfully provided on YAML nf pipeline github wiki). I then decompressed it and supplied the path to the directory
// The DB is very large (~1.9 GB) and I'm struggling how best to work around this.... apologies for the trouble caused to get this running!

include { METAPHLAN3_RUN } from '../../../../software/metaphlan3/run/main.nf' addParams( options: [ 'args':'--index mpa_v30_CHOCOPhlAn_201901' ] ) //database version (latest)

workflow test_metaphlan3_single_end {

    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    db    = channel.fromPath("/home/AD/mgordon/DATA/microbiome_testdata/testdata/YAMP_data/metaphlan_databases", type: 'dir', checkIfExists: true)

    METAPHLAN3_RUN ( input, db )
}

workflow test_metaphlan3_paired_end {

    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    db    = channel.fromPath("/home/AD/mgordon/DATA/microbiome_testdata/testdata/YAMP_data/metaphlan_databases", type: 'dir', checkIfExists: true)

    METAPHLAN3_RUN ( input, db )
}

workflow test_metaphlan3_sam {

    input = [ [ id:'test'], // meta map
              [ file("/home/AD/mgordon/PROJECTS/B034_NextFlow_Metagenomics/modules/tests/config/metaphlan3_testdata/toy_metagenome.sam", checkIfExists: true) ]
            ]

    db    = channel.fromPath("/home/AD/mgordon/DATA/microbiome_testdata/testdata/YAMP_data/metaphlan_databases", type: 'dir', checkIfExists: true)

    METAPHLAN3_RUN ( input, db )
}

workflow test_metaphlan3_fasta {

    input = [ [ id:'test'], // meta map
              [ file("/home/AD/mgordon/PROJECTS/B034_NextFlow_Metagenomics/modules/tests/config/metaphlan3_testdata/toy_metagenome.fasta.gz", checkIfExists: true) ]
            ]

    db    = channel.fromPath("/home/AD/mgordon/DATA/microbiome_testdata/testdata/YAMP_data/metaphlan_databases", type: 'dir', checkIfExists: true)

    METAPHLAN3_RUN ( input, db )
}


workflow test_metaphlan3_bowtie2out {

    input = [ [ id:'test'], // meta map
              [ file("/home/AD/mgordon/PROJECTS/B034_NextFlow_Metagenomics/modules/tests/config/metaphlan3_testdata/toy_metagenome.bowtie2out.txt", checkIfExists: true) ]
            ]

    db    = channel.fromPath("/home/AD/mgordon/DATA/microbiome_testdata/testdata/YAMP_data/metaphlan_databases", type: 'dir', checkIfExists: true )

    METAPHLAN3_RUN ( input, db )
}

