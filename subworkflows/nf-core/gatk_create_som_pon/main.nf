//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

params.mutect2_options      = [:]
params.gendbimport_options  = [:]
params.createsompon_options = [:]

include { GATK4_MUTECT2                     } from '../../../modules/gatk/mutect2/main'                     addParams( options: params.mutect2_options )
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/gatk/genomicsdbimport/main'            addParams( options: params.gendbimport_options )
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/gatk/createsomaticpanelofnormals/main' addParams( options: params.createsompon_options )

workflow GATK_CREATE_SOM_PON {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai/csi] ]

    path fasta
    path fastaidx
    path dict

    tuple val(meta) , path(input) , path(input_index) , []
    val  run_single
    val  run_pon
    val  run_mito

    tuple val(meta), path(vcf), path(tbi), path(intervalfile), [], []
    val run_intlist
    val run_updatewspace
    val input_map

    tuple val(meta), path(genomicsdb)

    main:
    ch_versions = Channel.empty()

    GATK4_MUTECT2 ( input , run_single , run_pon, run_mito , [] , fasta , fastaidx , dict , [], [] , [] , [] )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    GATK4_GENOMICSDBIMPORT GATK4_GENOMICSDBIMPORT ( input, run_intlist, run_updatewspace, input_map )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    GATK4_CREATESOMATICPANELOFNORMALS ( input, fasta, fastaidx, dict )
    ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions)

    emit:
    mutect2_vcf   = GATK4_MUTECT2.out.vcf                     // channel: [ val(meta), [ vcf ] ]
    mutect2_index = GATK4_MUTECT2.out.tbi                     // channel: [ val(meta), [ tbi ] ]
    mutect2_stats = GATK4_MUTECT2.out.stats                   // channel: [ val(meta), [ stats ] ]

    genomicsdb    = GATK4_GENOMICSDBIMPORT.out.genomicsdb     // channel: [ val(meta), [ genomicsdb ] ]

    pon_vcf       = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf // channel: [ val(meta), [ vcf.gz ] ]
    pon_index     = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi // channel: [ val(meta), [ tbi ] ]

    versions      = ch_versions                               // channel: [ versions.yml ]
}
