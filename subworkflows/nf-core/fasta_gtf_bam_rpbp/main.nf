//
// End-to-end Rp-Bp execution: build the Rp-Bp BED reference once, then run
// the per-sample ORF prediction chain on each Ribo-seq BAM. Deliberately
// avoids rpbp's `run-rpbp-pipeline` / `create-orf-profiles` wrappers - both
// invoke flexbar + bowtie + STAR on raw FASTQs, duplicating the upstream
// alignment supplied as input. Splitting the per-sample steps into separate
// processes also gives independent caching per step on resume.
//
// Order:
//   0. prepare-rpbp-genome (cohort-shared)        -> transcript / ORF BEDs
//   1. extract-metagene-profiles                  -> *.metagene-profile.csv.gz
//   2. estimate-metagene-profile-bayes-factors    -> *.metagene-periodicity-bayes-factors.csv.gz
//   3. select-periodic-offsets                    -> *.periodic-offsets.csv.gz
//   4. extract-orf-profiles                       -> *.profiles.mtx.gz
//   5. estimate-orf-bayes-factors                 -> *.bayes-factors.bed.gz
//   6. select-final-prediction-set                -> *.predicted-orfs.filtered.{bed.gz,dna.fa,protein.fa}
//

include { RPBP_PREPAREGENOME                } from '../../../modules/nf-core/rpbp/preparegenome/main'
include { RPBP_EXTRACTMETAGENEPROFILES      } from '../../../modules/nf-core/rpbp/extractmetageneprofiles/main'
include { RPBP_ESTIMATEMETAGENEBAYESFACTORS } from '../../../modules/nf-core/rpbp/estimatemetagenebayesfactors/main'
include { RPBP_SELECTPERIODICOFFSETS        } from '../../../modules/nf-core/rpbp/selectperiodicoffsets/main'
include { RPBP_EXTRACTORFPROFILES           } from '../../../modules/nf-core/rpbp/extractorfprofiles/main'
include { RPBP_ESTIMATEORFBAYESFACTORS      } from '../../../modules/nf-core/rpbp/estimateorfbayesfactors/main'
include { RPBP_SELECTFINALPREDICTIONSET     } from '../../../modules/nf-core/rpbp/selectfinalpredictionset/main'

workflow FASTA_GTF_BAM_RPBP {

    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ] - Ribo-seq BAMs
    ch_fasta_gtf  // channel (single value): [ val(meta), path(fasta), path(gtf) ]

    main:

    // 1. Prepare the Rp-Bp BED reference once per pipeline invocation.
    RPBP_PREPAREGENOME(ch_fasta_gtf)

    // 2. Take dedicated emissions so downstream sees them as proper files
    //    rather than reaching into the index dir (which loses the s3:// scheme
    //    under Fusion).
    ch_transcript_bed   = RPBP_PREPAREGENOME.out.transcript_bed  .map { _meta, bed -> bed }.first()
    ch_orfs_genomic_bed = RPBP_PREPAREGENOME.out.orfs_genomic_bed.map { _meta, bed -> bed }.first()
    ch_orfs_exons_bed   = RPBP_PREPAREGENOME.out.orfs_exons_bed  .map { _meta, bed -> bed }.first()
    ch_genome_fasta     = ch_fasta_gtf.map { _meta, fasta, _gtf -> fasta }.first()

    // 3. Per-sample ORF prediction chain.
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
    index            = RPBP_PREPAREGENOME.out.index                       // channel: [ val(meta), path(rpbp_index) ]
    transcript_bed   = RPBP_PREPAREGENOME.out.transcript_bed              // channel: [ val(meta), path(*.annotated.bed.gz) ]
    orfs_genomic_bed = RPBP_PREPAREGENOME.out.orfs_genomic_bed            // channel: [ val(meta), path(*.orfs-genomic.annotated.bed.gz) ]
    orfs_exons_bed   = RPBP_PREPAREGENOME.out.orfs_exons_bed              // channel: [ val(meta), path(*.orfs-exons.annotated.bed.gz) ]
    metagene         = RPBP_EXTRACTMETAGENEPROFILES.out.metagene          // channel: [ val(meta), path(*.metagene-profile.csv.gz) ]
    metagene_bf      = RPBP_ESTIMATEMETAGENEBAYESFACTORS.out.bayes_factors // channel: [ val(meta), path(*.metagene-periodicity-bayes-factors.csv.gz) ]
    periodic         = RPBP_SELECTPERIODICOFFSETS.out.periodic            // channel: [ val(meta), path(*.periodic-offsets.csv.gz) ]
    orf_profiles     = RPBP_EXTRACTORFPROFILES.out.profiles               // channel: [ val(meta), path(*.profiles.mtx.gz) ]
    lengths_offsets  = RPBP_EXTRACTORFPROFILES.out.lengths_offsets        // channel: [ val(meta), path(*.periodic_lengths_offsets.tsv) ]
    orf_bayes        = RPBP_ESTIMATEORFBAYESFACTORS.out.bayes_factors     // channel: [ val(meta), path(*.bayes-factors.bed.gz) ]
    predicted        = RPBP_SELECTFINALPREDICTIONSET.out.predicted        // channel: [ val(meta), path(*.predicted-orfs.filtered.bed.gz) ]
    dna_fasta        = RPBP_SELECTFINALPREDICTIONSET.out.dna_fasta        // channel: [ val(meta), path(*.predicted-orfs.filtered.dna.fa) ]
    protein_fasta    = RPBP_SELECTFINALPREDICTIONSET.out.protein_fasta    // channel: [ val(meta), path(*.predicted-orfs.filtered.protein.fa) ]
}
