//
// Rp-Bp post-alignment ORF prediction chain: one Nextflow process per rpbp tool.
//
// Replicates the chain that rpbp's `create-orf-profiles` + `predict-translated-orfs`
// wrappers would run internally, but as separate Nextflow processes so each step
// caches independently and only what changed is re-run. Avoids rpbp's internal
// flexbar/bowtie/STAR re-alignment - upstream alignment is supplied as the BAM.
//
// Order:
//   1. extract-metagene-profiles               -> *.metagene-profile.csv.gz
//   2. estimate-metagene-profile-bayes-factors -> *.metagene-periodicity-bayes-factors.csv.gz
//   3. select-periodic-offsets                 -> *.periodic-offsets.csv.gz
//   4. extract-orf-profiles                    -> *.profiles.mtx.gz
//   5. estimate-orf-bayes-factors              -> *.bayes-factors.bed.gz
//   6. select-final-prediction-set             -> *.predicted-orfs.filtered.{bed.gz,dna.fa,protein.fa}
//

include { RPBP_EXTRACTMETAGENEPROFILES      } from '../../../modules/nf-core/rpbp/extractmetageneprofiles/main'
include { RPBP_ESTIMATEMETAGENEBAYESFACTORS } from '../../../modules/nf-core/rpbp/estimatemetagenebayesfactors/main'
include { RPBP_SELECTPERIODICOFFSETS        } from '../../../modules/nf-core/rpbp/selectperiodicoffsets/main'
include { RPBP_EXTRACTORFPROFILES           } from '../../../modules/nf-core/rpbp/extractorfprofiles/main'
include { RPBP_ESTIMATEORFBAYESFACTORS      } from '../../../modules/nf-core/rpbp/estimateorfbayesfactors/main'
include { RPBP_SELECTFINALPREDICTIONSET     } from '../../../modules/nf-core/rpbp/selectfinalpredictionset/main'

workflow BAM_RPBP_PREDICTORFS {

    take:
    ch_bam               // channel: [ val(meta), path(bam), path(bai) ] - per Ribo-seq sample
    ch_transcript_bed    // channel: path - annotated transcripts BED from RPBP_PREPAREGENOME
    ch_orfs_genomic_bed  // channel: path - ORF genomic BED from RPBP_PREPAREGENOME
    ch_orfs_exons_bed    // channel: path - ORF exons BED from RPBP_PREPAREGENOME
    ch_genome_fasta      // channel: path - genome FASTA

    main:

    RPBP_EXTRACTMETAGENEPROFILES (
        ch_bam,
        ch_transcript_bed
    )

    RPBP_ESTIMATEMETAGENEBAYESFACTORS (
        RPBP_EXTRACTMETAGENEPROFILES.out.metagene
    )

    RPBP_SELECTPERIODICOFFSETS (
        RPBP_ESTIMATEMETAGENEBAYESFACTORS.out.bayes_factors
    )

    ch_extract_in = ch_bam
        .join(RPBP_SELECTPERIODICOFFSETS.out.periodic, by: 0)

    RPBP_EXTRACTORFPROFILES (
        ch_extract_in,
        ch_orfs_genomic_bed,
        ch_orfs_exons_bed
    )

    RPBP_ESTIMATEORFBAYESFACTORS (
        RPBP_EXTRACTORFPROFILES.out.profiles,
        ch_orfs_genomic_bed
    )

    RPBP_SELECTFINALPREDICTIONSET (
        RPBP_ESTIMATEORFBAYESFACTORS.out.bayes_factors,
        ch_genome_fasta
    )

    emit:
    metagene        = RPBP_EXTRACTMETAGENEPROFILES.out.metagene             // channel: [ val(meta), path(csv.gz) ]
    metagene_bf     = RPBP_ESTIMATEMETAGENEBAYESFACTORS.out.bayes_factors   // channel: [ val(meta), path(csv.gz) ]
    periodic        = RPBP_SELECTPERIODICOFFSETS.out.periodic               // channel: [ val(meta), path(csv.gz) ]
    orf_profiles    = RPBP_EXTRACTORFPROFILES.out.profiles                  // channel: [ val(meta), path(mtx.gz) ]
    lengths_offsets = RPBP_EXTRACTORFPROFILES.out.lengths_offsets           // channel: [ val(meta), path(tsv) ]
    orf_bayes       = RPBP_ESTIMATEORFBAYESFACTORS.out.bayes_factors        // channel: [ val(meta), path(bed.gz) ]
    predicted       = RPBP_SELECTFINALPREDICTIONSET.out.predicted           // channel: [ val(meta), path(bed.gz) ]
    dna_fasta       = RPBP_SELECTFINALPREDICTIONSET.out.dna_fasta           // channel: [ val(meta), path(fa) ]
    protein_fasta   = RPBP_SELECTFINALPREDICTIONSET.out.protein_fasta       // channel: [ val(meta), path(fa) ]
}
