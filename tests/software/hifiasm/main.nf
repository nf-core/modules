#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HIFIASM } from '../../../software/hifiasm/main.nf' addParams( options: [args:'-f0'] )

/* 
 * Test with long reads only
 */
workflow test_hifiasm_hifi_only {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/homo_sapiens/fastq/SRR10382244_mapped_to_contig.fastq", 
              checkIfExists: true) ] 
            ]
    hic = file('dummy_hic')
    paternal_kmer_dump = file('dummy_paternal_yak')
    maternal_kmer_dump = file('dummy_maternal_yak')
    HIFIASM ( input, hic, paternal_kmer_dump, maternal_kmer_dump, false, false)
}

/* 
 * Test with Hi-C reads for phasing
 */
// Hifiasm crashes with segfault on this test
// workflow test_hifiasm_with_hic {

//     def input = []
//     input = [ [ id:'test' ], // meta map
//               [ file("${launchDir}/tests/data/genomics/homo_sapiens/fastq/SRR10382244_mapped_to_contig.fastq", 
//               checkIfExists: true) ] 
//             ]
//     hic = Channel.fromFilePairs("${launchDir}/tests/data/genomics/homo_sapiens/fastq/SRR13061060_10000reads_mapped_{1,2}.fastq")
//     paternal_kmer_dump = file('dummy_paternal_yak')
//     maternal_kmer_dump = file('dummy_maternal_yak')
//     HIFIASM ( input, hic, paternal_kmer_dump, maternal_kmer_dump, true, false)
// }

/* 
 * Test with Hi-C reads for phasing
 */
workflow test_hifiasm_with_parental_reads {

    def input = []
    input = [ [ id:'test' ], // meta map
              [ file("${launchDir}/tests/data/genomics/homo_sapiens/fastq/SRR10382244_mapped_to_contig.fastq", 
              checkIfExists: true) ] 
            ]
    hic = file('dummy_hic')
    paternal_kmer_dump = file("${launchDir}/tests/data/genomics/homo_sapiens/yak/hg003_pat.yak")
    maternal_kmer_dump = file("${launchDir}/tests/data/genomics/homo_sapiens/yak/hg004_mat.yak")
    HIFIASM ( input, hic, paternal_kmer_dump, maternal_kmer_dump, false, true)
}
