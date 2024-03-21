include { METAPHLAN_MAKEDB      } from '../../../modules/nf-core/metaphlan/makedb/main'
include { METAPHLAN_METAPHLAN     } from '../../../modules/nf-core/metaphlan/metaphlan/main'
include { METAPHLAN_MERGEMETAPHLANTABLES     } from '../../../modules/nf-core/metaphlan/mergemetaphlantables/main'


workflow FASTQ_TAXONOMIC_PROFILE_METAPHLAN {

    take:
    ch_fastq

    main:

    ch_versions = Channel.empty()

    METAPHLAN_MAKEDB ( )
    ch_versions = ch_versions.mix(METAPHLAN_MAKEDB.out.versions.first())

    METAPHLAN_METAPHLAN ( ch_fastq, METAPHLAN_MAKEDB.out.db )
    ch_versions = ch_versions.mix(METAPHLAN_METAPHLAN.out.versions.first())

    METAPHLAN_MERGEMETAPHLANTABLES ( METAPHLAN_METAPHLAN.out.profile.map{ [ [id:'all_samples'], it[1] ] }.groupTuple( sort:true) )
    ch_versions = ch_versions.mix(METAPHLAN_MERGEMETAPHLANTABLES.out.versions.first())

    emit:
    merged_taxa      = METAPHLAN_MERGEMETAPHLANTABLES.out.txt


    versions = ch_versions
}

