//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

include { GATK4_MUTECT2                     } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createsomaticpanelofnormals/main'

workflow BAM_CREATE_SOM_PON_GATK {
    take:
    ch_mutect2_in       // channel: [ val(meta), path(input), path(input_index), path(interval_file) ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fai              // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]
    val_pon_norm        // string:  name for panel of normals
    ch_gendb_intervals  // channel: [ path(interval_file) ]

    main:
    ch_versions = Channel.empty()
    ch_input    = ch_mutect2_in

    //
    // Perform variant calling for each sample using mutect2 module in panel of normals mode.
    //
    GATK4_MUTECT2 (
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    //
    // Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    ch_vcf          = GATK4_MUTECT2.out.vcf.collect{it[1]}.toList()
    ch_index        = GATK4_MUTECT2.out.tbi.collect{it[1]}.toList()
    ch_dict_gendb   = ch_dict.map{meta, dict -> return dict}.toList()

    ch_gendb_input  = Channel.of([id:val_pon_norm])
        .combine(ch_vcf)
        .combine(ch_index)
        .combine(ch_gendb_intervals)
        .combine(ch_dict_gendb)
        .map{meta, vcf, tbi, interval, dict -> [meta, vcf, tbi, interval, [], dict]}

    GATK4_GENOMICSDBIMPORT ( ch_gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    //
    //Panel of normals made from genomicsdb workspace using createsomaticpanelofnormals.
    //
    GATK4_CREATESOMATICPANELOFNORMALS ( GATK4_GENOMICSDBIMPORT.out.genomicsdb, ch_fasta, ch_fai, ch_dict )
    ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions.first())

    emit:
    mutect2_vcf      = GATK4_MUTECT2.out.vcf                        // channel: [ val(meta), path(vcf) ]
    mutect2_index    = GATK4_MUTECT2.out.tbi                        // channel: [ val(meta), path(tbi) ]
    mutect2_stats    = GATK4_MUTECT2.out.stats                      // channel: [ val(meta), path(stats) ]
    genomicsdb       = GATK4_GENOMICSDBIMPORT.out.genomicsdb        // channel: [ val(meta), path(genomicsdb) ]
    pon_vcf          = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf    // channel: [ val(meta), path(vcf) ]
    pon_index        = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi    // channel: [ val(meta), path(tbi) ]

    versions         = ch_versions                                  // channel: [ path(versions.yml) ]
}
