//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//
params.mutect2_options      = [args: '--max-mnp-distance 0']
params.gendbimport_options  = [:]
params.createsompon_options = [:]

include { GATK4_MUTECT2                     } from '../../../modules/gatk4/mutect2/main'                     addParams( options: params.mutect2_options )
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/gatk4/genomicsdbimport/main'            addParams( options: params.gendbimport_options )
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/gatk4/createsomaticpanelofnormals/main' addParams( options: params.createsompon_options )

workflow GATK_CREATE_SOM_PON {
    take:
    ch_mutect2_in       // channel: [ val(meta), [ input ], [ input_index ], [] ]
    fasta               // channel: /path/to/reference/fasta
    fai                 // channel: /path/to/reference/fasta/index
    dict                // channel: /path/to/reference/fasta/dictionary
    pon_name            // channel: name for panel of normals
    interval_file       // channel: /path/to/interval/file

    main:
    ch_versions      = Channel.empty()
    input = channel.from(ch_mutect2_in)
    //
    //Perform variant calling for each sample using mutect2 module in panel of normals mode.
    //
    GATK4_MUTECT2 ( input, false, true, false, [], fasta, fai, dict, [], [], [], [] )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    //
    //Convert all sample vcfs into a genomicsdb workspace using genomicsdbimport.
    //
    ch_vcf = GATK4_MUTECT2.out.vcf.collect{it[1]}.toList()
    ch_index = GATK4_MUTECT2.out.tbi.collect{it[1]}.toList()
    gendb_input = Channel.of([[ id:pon_name ]]).combine(ch_vcf).combine(ch_index).combine([interval_file]).combine(['']).combine([dict])
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )
    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())

    //
    //Panel of normals made from genomicsdb workspace using createsomaticpanelofnormals.
    //
    GATK4_GENOMICSDBIMPORT.out.genomicsdb.view()
    GATK4_CREATESOMATICPANELOFNORMALS ( GATK4_GENOMICSDBIMPORT.out.genomicsdb, fasta, fai, dict )
    ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS.out.versions.first())

    emit:
    mutect2_vcf      = GATK4_MUTECT2.out.vcf.collect()                     // channel: [ val(meta), [ vcf ] ]
    mutect2_index    = GATK4_MUTECT2.out.tbi.collect()                     // channel: [ val(meta), [ tbi ] ]
    mutect2_stats    = GATK4_MUTECT2.out.stats.collect()                   // channel: [ val(meta), [ stats ] ]

    genomicsdb       = GATK4_GENOMICSDBIMPORT.out.genomicsdb               // channel: [ val(meta), [ genomicsdb ] ]

    pon_vcf          = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf           // channel: [ val(meta), [ vcf.gz ] ]
    pon_index        = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi           // channel: [ val(meta), [ tbi ] ]

    versions         = ch_versions                                         // channel: [ versions.yml ]
}
