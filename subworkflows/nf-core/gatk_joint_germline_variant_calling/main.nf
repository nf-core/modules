//
// Run GATK haplotypecaller for all input samples, merge them with genomicsdbimport, perform joint genotyping with genotypeGVCFS and recalibrate with variantrecalibrator & applyvqsr.
//

include { GATK4_HAPLOTYPECALLER     as HAPLOTYPECALLER     } from '../../../modules/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT    as GENOMICSDBIMPORT   } from '../../../modules/gatk4/genomicsdbimport/main'
// include { GATK4_GENOTYPEGVCFS       as GENOTYPEGVCFS       } from '../../../modules/gatk4/genotypegvcfs/main'
// include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR } from '../../../modules/gatk4/variantrecalibrator/main'
// include { GATK4_APPLYVQSR           as APPLYVQSR           } from '../../../modules/gatk4/applyvqsr/main'

workflow GATK_JOINT_GERMLINE_VARIANT_CALLING {
    take:
    // gatk_joint_germline_variant_calling
    // haplotypecaller
    input                     // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    sites
    sites_tbi
    interval_file
    // genomicsdbimport
    joint_id
    // genotypegvcfs
    // path  intervals_bed
    // variantrecalibrator
    // val allelespecific
    // tuple path(resvcfs), path(restbis), val(reslabels)
    // val annotation
    // val mode
    // val create_rscript
    // applyvqsr
    // val truthsensitivity

    main:
    ch_versions = Channel.empty()
    haplotc_input = channel.from(input)

    //
    //Perform variant calling using haplotypecaller module.
    //
    HAPLOTYPECALLER ( haplotc_input, fasta, fai, dict, sites, sites_tbi, interval_file )
    ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions.first())

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    ch_vcf = HAPLOTYPECALLER.out.vcf.collect{it[1]}.toList()
    ch_index = HAPLOTYPECALLER.out.tbi.collect{it[1]}.toList()
    additional_gendb_opts = [ interval_file, [], []]
    gendb_input = Channel.of([[ id:joint_id ]]).combine(ch_vcf).combine(ch_index).combine(additional_gendb_opts)
    GENOMICSDBIMPORT ( gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GENOMICSDBIMPORT.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    // ch_pileup = GETPILEUPSUMMARIES.out.table.collect()
    //[] is a placeholder for the optional input where the matched normal sample would be passed in for tumor-normal samples, which is not necessary for this workflow.
    // ch_pileup.add([])
    // CALCULATECONTAMINATION ( ch_pileup, true )
    // ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
    // tuple val(meta), path(gvcf), path(gvcf_index)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //
    // ch_vcf =           MUTECT2.out.vcf.collect()
    // ch_tbi =           MUTECT2.out.tbi.collect()
    // ch_stats =         MUTECT2.out.stats.collect()
    //[] is added as a placeholder for the optional input file artifact priors, which is only used for tumor-normal samples and therefor isn't needed in this workflow.
    // ch_stats.add([])
    // ch_segment =       CALCULATECONTAMINATION.out.segmentation.collect()
    // ch_contamination = CALCULATECONTAMINATION.out.contamination.collect()
    //[] is added as a placeholder for entering a contamination estimate value, which is not needed as this workflow uses the contamination table instead.
    // ch_contamination.add([])
    // ch_filtermutect_in = ch_vcf.combine(ch_tbi, by: 0).combine(ch_stats, by: 0).combine(ch_segment, by: 0).combine(ch_contamination, by: 0)
    // FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    // ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)
    // tuple val(meta), path(vcf) , path(tbi)

    // FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
    // ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)
    // tuple val(meta), path(vcf), path(tbi), path(recal), path(recalidx), path(tranches)

    emit:
    haplotc_vcf            = HAPLOTYPECALLER.out.vcf.collect()                             // channel: [ val(meta), [ vcf ] ]
    haplotc_index          = HAPLOTYPECALLER.out.tbi.collect()                             // channel: [ val(meta), [ tbi ] ]

    genomicsdb       = GENOMICSDBIMPORT.out.genomicsdb               // channel: [ val(meta), [ genomicsdb ] ]

    // contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), [ contamination ] ]
    // segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), [ segmentation ] ]

    // filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), [ vcf ] ]
    // filtered_index         = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), [ tbi ] ]
    // filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), [ stats ] ]

    versions               = ch_versions                                           // channel: [ versions.yml ]
}
