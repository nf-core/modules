#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD     } from '../../../../modules/nf-core/bowtie2/build/main.nf'
include { BWA_INDEX         } from '../../../../modules/nf-core/bwa/index/main.nf'
include { BWAMEM2_INDEX     } from '../../../../modules/nf-core/bwamem2/index/main.nf'
include { DRAGMAP_HASHTABLE } from '../../../../modules/nf-core/dragmap/hashtable/main.nf'
include { FASTQ_ALIGN_DNA   } from '../../../../subworkflows/nf-core/fastq_align_dna/main.nf'
include { SNAPALIGNER_INDEX } from '../../../../modules/nf-core/snapaligner/index/main.nf'

workflow test_fastq_align_bowtie2_SE {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BOWTIE2_BUILD ( fasta )
    FASTQ_ALIGN_DNA ( input, BOWTIE2_BUILD.out.index, "bowtie2", true )
}

workflow test_fastq_align_bowtie2_PE {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)

        ]
    ]
    fasta = [
        [ id:'test' ],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    BOWTIE2_BUILD ( fasta )
    FASTQ_ALIGN_DNA ( input, BOWTIE2_BUILD.out.index, "bowtie2", true )
}

workflow test_fastq_align_bwa_SE {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    BWA_INDEX ( [ [:], fasta ] )
    FASTQ_ALIGN_DNA ( input, BWA_INDEX.out.index, "bwamem", true )
}

workflow test_fastq_align_bwa_PE {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    BWA_INDEX ( [ [:], fasta ] )
    FASTQ_ALIGN_DNA ( input, BWA_INDEX.out.index, "bwamem", true )
}

workflow test_fastq_align_bwamem2_SE {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    BWAMEM2_INDEX ( [ [:], fasta ] )
    FASTQ_ALIGN_DNA ( input, BWAMEM2_INDEX.out.index, "bwamem2", true )
}

workflow test_fastq_align_bwamem2_PE {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    BWAMEM2_INDEX ( [ [:], fasta ] )
    FASTQ_ALIGN_DNA ( input, BWAMEM2_INDEX.out.index, "bwamem2", true )
}

workflow test_fastq_align_dragmap_SE {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    DRAGMAP_HASHTABLE ( fasta )
    FASTQ_ALIGN_DNA ( input, DRAGMAP_HASHTABLE.out.hashmap, "dragmap", true )
}

workflow test_fastq_align_dragmap_PE {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id:'test'],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]
    DRAGMAP_HASHTABLE ( fasta )
    FASTQ_ALIGN_DNA ( input, DRAGMAP_HASHTABLE.out.hashmap, "dragmap", true )
}

workflow test_fastq_align_snapaligner_SE {
    input = [
        [ id:'test', single_end:true ], // meta map
        [file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true)]
    ]
    fasta = [
        [id:"test"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        [],
        [],
        []
    ]
    SNAPALIGNER_INDEX ( fasta )
    FASTQ_ALIGN_DNA ( input, SNAPALIGNER_INDEX.out.index, "snap", true )
}

workflow test_fastq_align_snapaligner_PE {
    input = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]
    fasta = [
        [id:"test"],
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true),
        [],
        [],
        []
    ]
    SNAPALIGNER_INDEX ( fasta )
    FASTQ_ALIGN_DNA ( input, SNAPALIGNER_INDEX.out.index, "snap", true )
}
