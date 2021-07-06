#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../../software/star/genomegenerate/main.nf' addParams( options: [args: '--genomeSAindexNbases 9'] )
include { STAR_ALIGN          } from '../../../../software/star/align/main.nf'          addParams( options: [args: '--readFilesCommand zcat'] )

workflow test_star_alignment_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/generic/fastq/test_single_end.fastq.gz", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/generic/fasta/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf   = file("${launchDir}/tests/data/generic/gtf/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)
    
    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf )
}

workflow test_star_alignment_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/generic/fastq/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/generic/fastq/test_R2.fastq.gz", checkIfExists: true) ] 
            ]
    fasta = file("${launchDir}/tests/data/generic/fasta/GCF_000019425.1_ASM1942v1_genomic.fna", checkIfExists: true)
    gtf   = file("${launchDir}/tests/data/generic/gtf/GCF_000019425.1_ASM1942v1_genomic.gtf", checkIfExists: true)

    STAR_GENOMEGENERATE ( fasta, gtf )
    STAR_ALIGN ( input, STAR_GENOMEGENERATE.out.index, gtf )
}
