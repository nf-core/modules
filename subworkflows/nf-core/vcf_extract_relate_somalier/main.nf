include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_EXTRACT_RELATE_SOMALIER {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi), val(count) ]
        ch_fasta                // channel: [mandatory] [ path(fasta) ]
        ch_fasta_fai            // channel: [mandatory] [ path(fai) ]
        ch_somalier_sites       // channel: [mandatory] [ path(somalier_sites_vcf) ]
        ch_peds                 // channel: [mandatory] [ val(meta), path(ped) ]
        ch_sample_groups        // channel: [optional]  [ path(txt) ]
        val_common_id           // string:  [optional]  A common identifier for the samples that need to be related. Has to be given when using single sample VCFs
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

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_somalierextract_input = ch_input.no_tbi
        .join(TABIX_TABIX.out.tbi)
        .mix(ch_input.tbi)

    SOMALIER_EXTRACT(
        ch_somalierextract_input,
        ch_fasta,
        ch_fasta_fai,
        ch_somalier_sites
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions.first())

    ch_somalierrelate_input = SOMALIER_EXTRACT.out.extract
        .join(ch_vcfs, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, extract, vcf, tbi, count ->
            new_meta = val_common_id ? meta + [id:meta[val_common_id]] : meta
            [ count ? groupKey(new_meta, count): new_meta, extract ]
        }
        .groupTuple()
        .join(ch_peds, failOnDuplicate: true, failOnMismatch: true)
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
    extract        = SOMALIER_EXTRACT.out.extract       // channel: [ val(meta), path(extract) ]
    html           = SOMALIER_RELATE.out.html           // channel: [ val(meta), path(html) ]
    pairs_tsv      = SOMALIER_RELATE.out.pairs_tsv      // channel: [ val(meta), path(tsv) ]
    samples_tsv    = SOMALIER_RELATE.out.samples_tsv    // channel: [ val(meta), path(tsv) ]
    versions       = ch_versions                        // channel: [ path(versions.yml) ]
}
