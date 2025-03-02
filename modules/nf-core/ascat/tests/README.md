Pytest scripts added here for reference.

## Table of contents

- main.nf
- nextflow.config
- test.yml

## main.nf

```
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASCAT as ASCAT_SIMPLE }            from '../../../../modules/nf-core/ascat/main.nf'
include { ASCAT as ASCAT_PLOIDY_AND_PURITY } from '../../../../modules/nf-core/ascat/main.nf'
include { ASCAT as ASCAT_CRAM }              from '../../../../modules/nf-core/ascat/main.nf'
include { UNZIP as UNZIP_ALLELES }           from '../../../../modules/nf-core/unzip/main.nf'
include { UNZIP as UNZIP_LOCI }              from '../../../../modules/nf-core/unzip/main.nf'
include { UNZIP as UNZIP_GC }                from '../../../../modules/nf-core/unzip/main.nf'
include { UNZIP as UNZIP_RT }                from '../../../../modules/nf-core/unzip/main.nf'

workflow test_ascat {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    ASCAT_SIMPLE ( input , [], [], [], [], [], [])
}

// extended tests running with 1000 genomes data. Data is downloaded as follows:
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00154/alignment/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00154/alignment/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00155/alignment/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00155/alignment/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai
wget https://zenodo.org/records/14008443/files/G1000_loci_WGS_hg38.zip
unzip G1000_loci_WGS_hg38.zip
wget https://zenodo.org/records/14008443/files/G1000_alleles_WGS_hg38.zip
unzip G1000_alleles_WGS_hg38.zip
// wget https://www.dropbox.com/s/v0tgr1esyoh1krw/GC_G1000_hg19.zip
// wget https://www.dropbox.com/s/50n7xb06x318tgl/RT_G1000_hg19.zip


samtools view test2.paired_end.markduplicates.sorted.bam | awk '{if (NR==1) {min=$4; max=$4} if ($4<min) {min=$4} if ($4>max) {max=$4}} END {print "Min position:", min, "\nMax position:", max}'


// workflow test_ascat_with_ploidy_and_purity {
//     input = [ [ id:'test', single_end:false ], // meta map
//             file("HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
//             file("HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true),
//             file("HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
//             file("HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true)
//     ]

//     allele_path  = file("G1000_alleles_hg19.zip", checkIfExists: true)
//     allele_files = [[ id: allele_path.BaseName ], allele_path  ]

//     loci_path    = file("G1000_loci_hg19.zip", checkIfExists: true)
//     loci_files   = [[ id: loci_path.BaseName ], loci_path  ]

//     gc_path      = file("GC_G1000_hg19.zip", checkIfExists: true)
//     gc_file      = [[ id: gc_path.BaseName ], gc_path  ]

//     rt_path     = file("RT_G1000_hg19.zip", checkIfExists: true)
//     rt_file     = [[ id: rt_path.BaseName ], rt_path  ]

//     UNZIP_ALLELES(allele_files)
//     UNZIP_LOCI(loci_files)
//     UNZIP_GC(gc_file)

//     ASCAT_PLOIDY_AND_PURITY ( input ,
//                                 UNZIP_ALLELES.out.unzipped_archive.map{ it[1] },
//                                 UNZIP_LOCI.out.unzipped_archive.map{ it[1] },
//                                 [],                                             // optional bed_file for WES
//                                 [],                                             // optional fasta
//                                 UNZIP_GC.out.unzipped_archive.map{ it[1] },     // optional GC_correction
//                                 [])                                             // optional RT_correction
// }

// extended tests running with 1000 genomes data. Data is downloaded as follows:
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00145/alignment/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00145/alignment/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00146/alignment/HG00146.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00146/alignment/HG00146.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram
// wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
// wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

// workflow test_ascat_with_crams {
//     input = [
//        [ id:'test', single_end:false ], // meta map
//         file("HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram", checkIfExists: true),
//         file("HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai", checkIfExists: true),
//         file("HG00146.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram", checkIfExists: true),
//         file("HG00146.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai", checkIfExists: true)
//     ]

//     allele_path  = file("G1000_alleles_hg19.zip", checkIfExists: true)
//     allele_files = [[ id: allele_path.BaseName ], allele_path  ]

//     loci_path    = file("G1000_loci_hg19.zip", checkIfExists: true)
//     loci_files   = [[ id: loci_path.BaseName ], loci_path  ]

//     gc_path      = file("GC_G1000_hg19.zip", checkIfExists: true)
//     gc_file      = [[ id: gc_path.BaseName ], gc_path  ]

//     rt_path     = file("RT_G1000_hg19.zip", checkIfExists: true)
//     rt_file     = [[ id: rt_path.BaseName ], rt_path  ]

//     fasta   = file("human_g1k_v37.fasta", checkIfExists: true)

//     UNZIP_ALLELES(allele_files)
//     UNZIP_LOCI(loci_files)
//     UNZIP_GC(gc_file)
//     UNZIP_RT(rt_file)

//     ASCAT_CRAM ( input ,
//                 UNZIP_ALLELES.out.unzipped_archive.map{ it[1] },
//                 UNZIP_LOCI.out.unzipped_archive.map{ it[1] },
//                 [],
//                 fasta,
//                 UNZIP_GC.out.unzipped_archive.map{ it[1] },
//                 UNZIP_RT.out.unzipped_archive.map{ it[1] })
// }
```

## nextflow.config

```
process {

    publishDir = { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" }


    withName: ASCAT_SIMPLE {
        ext.args = [
          gender    :  'XY',
          genomeVersion : 'hg19',
          minCounts :  '1',
          min_base_qual : '1',
          min_map_qual  : '1',
          chrom_names  : 'c("21","22")'
          ]
    }



    withName: ASCAT_PLOIDY_AND_PURITY {
        ext.args = [
          gender : 'XX',
          genomeVersion : 'hg19',
          ploidy  : '1.7',
          purity : '0.24',
          chrom_names  : 'c("21","22")',

        ]
    }

    withName: ASCAT_CRAM {
        ext.args = [
          gender    : 'XX',
          genomeVersion : 'hg19',
          ref_fasta : '/mnt/volume/ascat/human_g1k_v37.fasta',
          chrom_names  : 'c("21","22")'
        ]
    }

}
```

## test.yml

```
- name: ascat test_ascat
  command: nextflow run ./tests/modules/nf-core/ascat -entry test_ascat -c ./tests/config/nextflow.config -stub-run
  tags:
    - ascat
  files:
    - path: output/ascat/test.after_correction.gc_rt.test.tumour.germline.png
    - path: output/ascat/test.after_correction.gc_rt.test.tumour.tumour.png
    - path: output/ascat/test.before_correction.test.tumour.germline.png
    - path: output/ascat/test.before_correction.test.tumour.tumour.png
    - path: output/ascat/test.cnvs.txt
    - path: output/ascat/test.metrics.txt
    - path: output/ascat/test.normal_alleleFrequencies_chr21.txt
    - path: output/ascat/test.normal_alleleFrequencies_chr22.txt
    - path: output/ascat/test.purityploidy.txt
    - path: output/ascat/test.segments.txt
    - path: output/ascat/test.tumour.ASPCF.png
    - path: output/ascat/test.tumour.sunrise.png
    - path: output/ascat/test.tumour_alleleFrequencies_chr21.txt
    - path: output/ascat/test.tumour_alleleFrequencies_chr22.txt
    - path: output/ascat/test.tumour_normalBAF.txt
    - path: output/ascat/test.tumour_normalLogR.txt
    - path: output/ascat/test.tumour_tumourBAF.txt
    - path: output/ascat/test.tumour_tumourLogR.txt
# - name: ascat test_ascat_with_ploidy_and_purity
#   command: nextflow run ./tests/modules/nf-core/ascat -entry test_ascat_with_ploidy_and_purity -c ./tests/config/nextflow.config
#   tags:
#     - ascat
#     - bam
#   files:
#   - path: output/ascat/test.normal_alleleFrequencies_chr22.txt
#     md5sum: 34a7ff0d88d232e6dca236f74bd616a3
#   - path: output/ascat/test.tumour_normalLogR.txt
#     md5sum: ad3583bb9f02ce371b3c6e738e01253c
#   - path: output/ascat/test.segments.txt
#     md5sum: 0eba21b2c35bd6b20422b5d3f0d5268a
#   - path: output/ascat/test.tumour_normalBAF.txt
#     md5sum: 0d6ce73d3220237681bf24ec4119e1fa
#   - path: output/ascat/test.tumour_tumourBAF.txt
#     md5sum: 80eaf8b11b072523f8cc5f1c0810df09
#   - path: output/ascat/test.after_correction_gc.test.tumour.germline.png
#     md5sum: a26376403f8118a391f0408d9260b8cd
#   - path: output/ascat/test.cnvs.txt
#     md5sum: 4f99f706b2e498182bbddfaf8b0fafda
#   - path: output/ascat/test.normal_alleleFrequencies_chr21.txt
#     md5sum: 83665d6887709676ed28a3df9398cfed
#   - path: output/ascat/test.tumour.sunrise.png
#     md5sum: 025d911f77cc48e02be88c027a0f5c13
#   - path: output/ascat/test.tumour_alleleFrequencies_chr21.txt
#     md5sum: 3774b4cb495816ee481b907c3e515f19
#   - path: output/ascat/test.before_correction.test.tumour.germline.png
#     md5sum: a26376403f8118a391f0408d9260b8cd
#   - path: output/ascat/test.metrics.txt
#     md5sum: 0983a53cc0db5918b55d04d9e8a3d96a
#   - path: output/ascat/test.tumour.ASCATprofile.png
#     md5sum: 7a90b51e32ee779f18979217416475b9
#   - path: output/ascat/test.purityploidy.txt
#     md5sum: 48fa643c43bf81e76668f58c0fcf4cdd
#   - path: output/ascat/test.before_correction.test.tumour.tumour.png
#     md5sum: 30b9ead13a68cfa6a68d1e8512383c64
#   - path: output/ascat/test.tumour.ASPCF.png
#     md5sum: 071538088f2e41d335aa530761c5e1e9
#   - path: output/ascat/test.tumour.rawprofile.png
#     md5sum: 64e700602181e5323c1ff6d5cc48b29a
#   - path: output/ascat/test.after_correction_gc.test.tumour.tumour.png
#     md5sum: f3f7bcda5023092a0facb826099f134e
#   - path: output/ascat/test.tumour_tumourLogR.txt
#     md5sum: 1a7add5dd2697dbee1ec83c887ff0a81
#   - path: output/ascat/test.tumour_alleleFrequencies_chr22.txt
#     md5sum: d8e1eed4d6d1a91532adac1ce5a111cb

# - name: ascat test_ascat_with_ploidy_and_purity
#   command: nextflow run ./tests/modules/nf-core/ascat -entry test_ascat_with_crams -c ./tests/config/nextflow.config
#   tags:
#     - ascat
#     - cram
#   files:
# - path: output/ascat/test.after_correction_gc_rt.test.tumour.tumour.png
#   md5sum: 04fb5a06abaa64eeb6008f3ff321b806
# - path: output/ascat/test.tumour_alleleFrequencies_chr22.txt
#   md5sum: 17058c76cc3772363700ce6b604ad4e9
# - path: output/ascat/test.tumour.sunrise.png
#   md5sum: 1c933eb599628c7ae6fc64cf2093f2fb
# - path: output/ascat/test.after_correction_gc_rt.test.tumour.germline.png
#   md5sum: 1e29c586cff82140c653232da5f780e6
# - path: output/ascat/test.before_correction.test.tumour.germline.png
#   md5sum: 1e29c586cff82140c653232da5f780e6
# - path: output/ascat/test.tumour_normalBAF.txt
#   md5sum: 2e5c32f9ba7a573974c55abb59ae7c7d
# - path: output/ascat/test.normal_alleleFrequencies_chr21.txt
#   md5sum: 349d516eb745002265cee66fe9ff3c72
# - path: output/ascat/test.tumour_normalLogR.txt
#   md5sum: 557b158d19f1b3c678a7e2638c1a595a
# - path: output/ascat/test.metrics.txt
#   md5sum: 59b6db8eea49d45845d03962161e1747
# - path: output/ascat/test.tumour_alleleFrequencies_chr21.txt
#   md5sum: 60637b7eba083c3f9af3e37055bae43b
# - path: output/ascat/test.cnvs.txt
#   md5sum: 68b329da9893e34099c7d8ad5cb9c940
# - path: output/ascat/test.segments.txt
#   md5sum: 68b329da9893e34099c7d8ad5cb9c940
# - path: output/ascat/test.before_correction.test.tumour.tumour.png
#   md5sum: 81397df24928ff2aa6e6b99475f2f580
# - path: output/ascat/test.tumour_tumourBAF.txt
#   md5sum: c63500836886125aa0450691377d857e
# - path: output/ascat/test.tumour.ASPCF.png
#   md5sum: d2ca81a06d65fbe491894612c9ebc238
# - path: output/ascat/test.normal_alleleFrequencies_chr22.txt
#   md5sum: e439014caab3a0f03efb563741d5da08
# - path: output/ascat/test.tumour_tumourLogR.txt
#   md5sum: ebcd14ecdaaecf47e0d647d3871d3739
# - path: output/ascat/test.purityploidy.txt
#   md5sum: f1484c2b120834d3db8774ad02a038b
```
