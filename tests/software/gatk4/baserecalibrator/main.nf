#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BASERECALIBRATOR } from '../../../../software/gatk4/baserecalibrator/main.nf' addParams( options: [:] )

workflow test_gatk4_baserecalibrator {

  def input = []
  input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
  fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
  fastaidx = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
  dict = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.dict", checkIfExists: true)

  sites = file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz", checkIfExists: true)
  sites_tbi = file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz.tbi", checkIfExists: true)
  GATK4_BASERECALIBRATOR ( input, fasta, fastaidx, dict, [], sites, sites_tbi )
}

workflow test_gatk4_baserecalibrator_intervals {

  def input = []
  input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
  fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
  fastaidx = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
  dict = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.dict", checkIfExists: true)
  intervals = file("${launchDir}/tests/data/genomics/sarscov2/bed/test.bed", checkIfExists: true)
  sites = file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz", checkIfExists: true)
  sites_tbi = file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz.tbi", checkIfExists: true)
  GATK4_BASERECALIBRATOR ( input, fasta, fastaidx, dict, intervals, sites, sites_tbi )
}


workflow test_gatk4_baserecalibrator_multiple_sites {

  def input = []
  input = [ [ id:'test' ], // meta map
            file("${launchDir}/tests/data/genomics/sarscov2/bam/test_paired_end.sorted.bam", checkIfExists: true) ]
  fasta = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta", checkIfExists: true)
  fastaidx = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.fasta.fai", checkIfExists: true)
  dict = file("${launchDir}/tests/data/genomics/sarscov2/fasta/test_genome.dict", checkIfExists: true)

  sites = [file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz", checkIfExists: true),
           file("${launchDir}/tests/data/genomics/sarscov2/vcf/test2.vcf.gz", checkIfExists: true)]
  sites_tbi = [file("${launchDir}/tests/data/genomics/sarscov2/vcf/test.vcf.gz.tbi", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/sarscov2/vcf/test2.vcf.gz.tbi", checkIfExists: true)]
  GATK4_BASERECALIBRATOR ( input, fasta, fastaidx, dict, [], sites, sites_tbi )
}
