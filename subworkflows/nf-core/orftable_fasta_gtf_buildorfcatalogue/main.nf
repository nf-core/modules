//
// Build a multi-caller Ribo-seq ORF catalogue from per-sample, per-caller
// ORF prediction tables. Composes the normaliser (run once per caller-tagged
// input), the merger (one cohort-level invocation), and bedtools/getfasta +
// seqkit/translate to produce the catalogue AA FASTA.
//

include { CUSTOM_ORFNORMALISE } from '../../../modules/nf-core/custom/orfnormalise/main'
include { CUSTOM_ORFMERGE     } from '../../../modules/nf-core/custom/orfmerge/main'
include { BEDTOOLS_GETFASTA   } from '../../../modules/nf-core/bedtools/getfasta/main'
include { SEQKIT_TRANSLATE    } from '../../../modules/nf-core/seqkit/translate/main'

workflow ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE {

    take:
    ch_orf_tables  // channel: [ val(meta), path(orf_table), val(caller) ]
                   //          caller in {ribocode, ribotish, ribotricer, rpbp, price};
                   //          all caller outputs flow through one channel with the
                   //          caller id carried as a per-record val (not in meta).
    ch_fasta       // channel: [ val(meta), path(fasta) ]   - reference genome FASTA
    ch_gtf         // channel: [ val(meta), path(gtf)   ]   - reference GTF (used by
                   //          ribocode/ribotish normalisers; ignored by rpbp/price)

    main:

    // 1. Normalise each caller's output. The same module is invoked once per
    //    input emission; dispatch happens inside the template based on the
    //    `caller` val. Append the caller to meta.id so the normaliser's
    //    default `${meta.id}` prefix yields caller-disambiguated filenames
    //    (the merger stages multiple BED12s into `beds/*` and needs unique
    //    names per caller-sample combination).
    ch_normalise_in = ch_orf_tables.map { meta, table, caller ->
        [ meta + [ id: "${meta.id}.${caller}" ], table, caller ]
    }
    CUSTOM_ORFNORMALISE ( ch_normalise_in, ch_gtf.first() )

    // 2. Gather all normalised BED12s + sidecar TSVs across callers and
    //    samples into a single cohort-keyed channel. `.collect()` on an
    //    empty upstream channel still emits `[]`, so filter out that case
    //    rather than feeding the merger a zero-arity input list.
    ch_all_beds = CUSTOM_ORFNORMALISE.out.bed12
        .map { _meta, bed -> bed }
        .collect()
        .filter { it.size() > 0 }
        .map { beds -> [ 'cohort', beds ] }

    ch_all_tsvs = CUSTOM_ORFNORMALISE.out.tsv
        .map { _meta, tsv -> tsv }
        .collect()
        .filter { it.size() > 0 }
        .map { tsvs -> [ 'cohort', tsvs ] }

    ch_merge_in = ch_all_beds
        .combine(ch_all_tsvs, by: 0)
        .map { _key, beds, tsvs -> [ [ id: 'cohort' ], beds, tsvs ] }

    CUSTOM_ORFMERGE ( ch_merge_in )

    // 3. Lift the merged BED12 into nucleotide then amino-acid FASTA.
    //    `bedtools getfasta -split` walks BED12 blocks in mRNA order;
    //    `seqkit translate --trim` drops trailing stops.
    BEDTOOLS_GETFASTA ( CUSTOM_ORFMERGE.out.bed12, ch_fasta.map { _meta, fa -> fa }.first() )
    SEQKIT_TRANSLATE  ( BEDTOOLS_GETFASTA.out.fasta )

    emit:
    catalogue_bed12    = CUSTOM_ORFMERGE.out.bed12
    catalogue_tsv      = CUSTOM_ORFMERGE.out.catalogue_tsv
    orf_to_gene_tsv    = CUSTOM_ORFMERGE.out.orf_to_gene_tsv
    catalogue_aa_fasta = SEQKIT_TRANSLATE.out.fastx
    multiqc            = CUSTOM_ORFMERGE.out.multiqc
}
