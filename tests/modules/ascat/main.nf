#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASCAT as ASCAT_SIMPLE }            from '../../../modules/ascat/main.nf'
include { ASCAT as ASCAT_PLOIDY_AND_PURITY } from '../../../modules/ascat/main.nf'
include { ASCAT as ASCAT_CRAM }              from '../../../modules/ascat/main.nf'
include { UNZIP as UNZIP_ALLELES }           from '../../../modules/unzip/main.nf'
include { UNZIP as UNZIP_LOCI }              from '../../../modules/unzip/main.nf'




workflow test_ascat {
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    ASCAT_SIMPLE ( input , [], [])
}





// extended tests running with 1000 genomes data. Data is downloaded as follows:
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00154/alignment/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00154/alignment/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00155/alignment/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam
// wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00155/alignment/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai
workflow test_ascat_with_ploidy_and_purity {  
    input = [ [ id:'test', single_end:false ], // meta map
         file("/mnt/volume/ascat/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
         file("/mnt/volume/ascat/HG00154.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true),
         file("/mnt/volume/ascat/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam", checkIfExists: true),
         file("/mnt/volume/ascat/HG00155.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam.bai", checkIfExists: true)
    ]

    allele_path  = file("/mnt/volume/repos/modules/test_ascat2/G1000_alleles_hg19.zip", checkIfExists: true)   
    allele_files = [[ id: allele_path.BaseName ], allele_path  ]

    loci_path    = file("/mnt/volume/repos/modules/test_ascat2/G1000_loci_hg19.zip", checkIfExists: true)   
    loci_files   = [[ id: loci_path.BaseName ], loci_path  ]

    // nextflow run ../tests/modules/ascat -entry test_ascat_with_ploidy_and_purity -c ../tests/config/nextflow.config ../tests/ascat/nextflow.config
    //loci_files   = "https://www.dropbox.com/s/l3m0yvyca86lpwb/G1000_loci_hg19.zip"
    //"/mnt/volume/ascat/alleles/G1000_alleles_hg19_chr"
    //"/mnt/volume/ascat/loci/G1000_loci_hg19_chr"

    UNZIP_ALLELES(allele_files)
    UNZIP_LOCI(loci_files)
    //UNZIP_ALLELES.out.unzipped_archive.view()
    //UNZIP_LOCI.out.unzipped_archive.view()
    ASCAT_PLOIDY_AND_PURITY ( input , UNZIP_ALLELES.out.unzipped_archive.map{ it[1] } , UNZIP_LOCI.out.unzipped_archive.map{ it[1] } )

}

// extended tests running with 1000 genomes data. Data is downloaded as follows:
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00145/alignment/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00145/alignment/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00146/alignment/HG00146.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai
// wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00146/alignment/HG00146.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram
//workflow test_ascat_with_crams {
//    input = [
//        [ id:'test', single_end:false ], // meta map
//        file("/home/ec2-user/input_files/crams/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram", checkIfExists: true),
//        file("/home/ec2-user/input_files/crams/HG00145.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram.crai", checkIfExists: true),
//        file("/home/ec2-user/input_files/crams/duplicate_test.cram", checkIfExists: true),
//        file("/home/ec2-user/input_files/crams/duplicate_test.cram.crai", checkIfExists: true)
//    ]
//
//    ASCAT_CRAM ( input , "/home/ec2-user/input_files/allele_files/G1000_alleles_hg19_chr", "/home/ec2-user/input_files/loci_files/G1000_alleles_hg19_chr")
//}



