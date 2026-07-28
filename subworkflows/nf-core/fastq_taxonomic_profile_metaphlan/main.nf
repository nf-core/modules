include { METAPHLAN_MAKEDB               } from '../../../modules/nf-core/metaphlan/makedb/main'
include { METAPHLAN_METAPHLAN            } from '../../../modules/nf-core/metaphlan/metaphlan/main'
include { METAPHLAN_MERGEMETAPHLANTABLES } from '../../../modules/nf-core/metaphlan/mergemetaphlantables/main'


workflow FASTQ_TAXONOMIC_PROFILE_METAPHLAN {
    take:
    ch_fastq

    main:

    METAPHLAN_MAKEDB()

    METAPHLAN_METAPHLAN(ch_fastq, METAPHLAN_MAKEDB.out.db, false)

    metaphlan_merged_profiles_txt = METAPHLAN_MERGEMETAPHLANTABLES(METAPHLAN_METAPHLAN.out.profile
        .map { _meta, profile -> [[id: 'all_samples'], profile] }
        .groupTuple(sort: { profile -> profile.getName() }))
        .txt

    emit:
    merged_taxa = metaphlan_merged_profiles_txt
}
