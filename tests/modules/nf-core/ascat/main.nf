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
// wget https://www.dropbox.com/s/l3m0yvyca86lpwb/G1000_loci_hg19.zip
// wget https://www.dropbox.com/s/3fzvir3uqe3073d/G1000_alleles_hg19.zip
// wget https://www.dropbox.com/s/v0tgr1esyoh1krw/GC_G1000_hg19.zip
// wget https://www.dropbox.com/s/50n7xb06x318tgl/RT_G1000_hg19.zip

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



