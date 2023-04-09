include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_EXTRACT_RELATE_SOMALIER {
    take:
        ch_vcfs                 // channel: [mandatory] [ meta, vcf, tbi, count ]
        ch_fasta                // channel: [mandatory] [ fasta ]
        ch_fasta_fai            // channel: [mandatory] [ fai ]
        ch_somalier_sites       // channel: [mandatory] [ somalier_sites_vcf ]
        ch_peds                 // channel: [mandatory] [ meta, ped ]
        ch_sample_groups        // channel: [optional]  [ txt ]
    main:

    ch_versions         = Channel.empty()

    ch_input = ch_vcfs
        .branch { meta, vcf, tbi, count ->
            tbi: tbi != []
                return [ meta, vcf, tbi ]
            no_tbi: tbi == []
                return [ meta, vcf ]
        }

    TABIX_TABIX(
        ch_input.no_tbi
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    ch_somalierextract_input = ch_input.no_tbi
        .join(TABIX_TABIX.out.tbi)
        .mix(ch_input.tbi)

    SOMALIER_EXTRACT(
        ch_somalierextract_input,
        ch_fasta,
        ch_fasta_fai,
        ch_somalier_sites
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions)

    ch_somalierrelate_input = SOMALIER_EXTRACT.out.extract
        .join(ch_vcfs)
        .map { meta, extract, vcf, tbi, count ->
            [ count ? groupKey(meta, count): meta, extract ]
        }
        .groupTuple()
        .join(ch_peds)
        .map { meta, extract, ped ->
            extract2 = extract[0] instanceof ArrayList ? extract[0] : extract
            [ meta, extract2, ped ]
        }

    SOMALIER_RELATE(
        ch_somalierrelate_input,
        ch_sample_groups
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions)

    emit:
    extract        = SOMALIER_EXTRACT.out.extract
    html           = SOMALIER_RELATE.out.html
    pairs_tsv      = SOMALIER_RELATE.out.pairs_tsv
    samples_tsv    = SOMALIER_RELATE.out.samples_tsv
    versions       = ch_versions
}
