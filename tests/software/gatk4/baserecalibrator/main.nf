#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BASERECALIBRATOR } from '../../../../software/gatk4/baserecalibrator/main.nf' addParams( options: [:] )

workflow test_gatk4_baserecalibrator {
  input     = [ [ id:'test' ], // meta map
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
              ]
  fasta     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
  fai       = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
  dict      = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)
  sites     = file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz", checkIfExists: true)
  sites_tbi = file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz.tbi", checkIfExists: true)

  GATK4_BASERECALIBRATOR ( input, fasta, fai, dict, [], sites, sites_tbi )
}

workflow test_gatk4_baserecalibrator_intervals {
  input =   [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
            ]
  fasta     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
  fai       = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
  dict      = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)
  intervals = file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true)
  sites     = file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz", checkIfExists: true)
  sites_tbi = file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz.tbi", checkIfExists: true)

  GATK4_BASERECALIBRATOR ( input, fasta, fai, dict, intervals, sites, sites_tbi )
}

workflow test_gatk4_baserecalibrator_multiple_sites {
  input     = [ [ id:'test' ], // meta map
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/bam/test_paired_end.sorted.bam", checkIfExists: true) 
              ]
  fasta     = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta", checkIfExists: true)
  fai       = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.fasta.fai", checkIfExists: true)
  dict      = file("${launchDir}/tests/data/genomics/sarscov2/genome/genome.dict", checkIfExists: true)
  sites     = [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test2.vcf.gz", checkIfExists: true)
              ]
  sites_tbi = [ file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz.tbi", checkIfExists: true),
                file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test2.vcf.gz.tbi", checkIfExists: true)
              ]

  GATK4_BASERECALIBRATOR ( input, fasta, fai, dict, [], sites, sites_tbi )
}
