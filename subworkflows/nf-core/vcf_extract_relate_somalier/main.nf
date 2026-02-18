include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'
include { TABIX_TABIX      } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_EXTRACT_RELATE_SOMALIER {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi), val(count) ]
        ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fasta_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_somalier_sites       // channel: [mandatory] [ path(somalier_sites_vcf) ]
        ch_peds                 // channel: [mandatory] [ val(meta), path(ped) ]
        ch_sample_groups        // channel: [optional]  [ path(txt) ]
        val_common_id           // string:  [optional]  A common identifier for the samples that need to be related. Has to be given when using single sample VCFs
    main:

    ch_input = ch_vcfs
        .branch { meta, vcf, tbi, _count ->
            tbi: tbi != []
                return [ meta, vcf, tbi ]
            no_tbi: tbi == []
                return [ meta, vcf ]
        }

    TABIX_TABIX(
        ch_input.no_tbi
    )

    ch_somalierextract_input = ch_input.no_tbi
        .join(TABIX_TABIX.out.index)
        .mix(ch_input.tbi)

    SOMALIER_EXTRACT(
        ch_somalierextract_input,
        ch_fasta,
        ch_fasta_fai,
        ch_somalier_sites
    )

    ch_somalierrelate_input = SOMALIER_EXTRACT.out.extract
        .join(ch_vcfs, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, extract, _vcf, _tbi, count ->
            def new_meta = val_common_id ? meta + [id:meta[val_common_id]] : meta
            [ count ? groupKey(new_meta, count): new_meta, extract ]
        }
        .groupTuple()
        .join(ch_peds, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, extract, ped ->
            def extract2 = extract[0] instanceof ArrayList ? extract[0] : extract
            def sorted_extract = extract2.sort { a, b -> file(a).name <=> file(b).name }
            // Check if meta is a GroupKey by checking for 'target' property
            def new_meta = meta.hasProperty('target') ? meta.target : meta
            [ new_meta, sorted_extract, ped ]
        } // Sort and flatten the extract list, remove the GroupKey wrapper if present

    SOMALIER_RELATE(
        ch_somalierrelate_input,
        ch_sample_groups
    )

    emit:
    extract        = SOMALIER_EXTRACT.out.extract       // channel: [ val(meta), path(extract) ]
    html           = SOMALIER_RELATE.out.html           // channel: [ val(meta), path(html) ]
    pairs_tsv      = SOMALIER_RELATE.out.pairs_tsv      // channel: [ val(meta), path(tsv) ]
    samples_tsv    = SOMALIER_RELATE.out.samples_tsv    // channel: [ val(meta), path(tsv) ]
}
