//
// End-to-end Rp-Bp execution: build the Rp-Bp BED reference once, then run
// the per-sample ORF prediction chain on each Ribo-seq BAM. Deliberately
// avoids rpbp's `run-rpbp-pipeline` / `create-orf-profiles` wrappers - both
// invoke flexbar + bowtie + STAR on raw FASTQs, duplicating the upstream
// alignment supplied as input. Splitting also gives independent caching
// per step on resume.
//

include { RPBP_PREPAREGENOME     } from '../../../modules/nf-core/rpbp/preparegenome/main'
include { BAM_RPBP_PREDICTORFS   } from '../bam_rpbp_predictorfs/main'

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
    BAM_RPBP_PREDICTORFS (
        ch_bam,
        ch_transcript_bed,
        ch_orfs_genomic_bed,
        ch_orfs_exons_bed,
        ch_genome_fasta
    )

    emit:
    index            = RPBP_PREPAREGENOME.out.index             // channel: [ val(meta), path(rpbp_index) ]
    transcript_bed   = RPBP_PREPAREGENOME.out.transcript_bed    // channel: [ val(meta), path(*.annotated.bed.gz) ]
    orfs_genomic_bed = RPBP_PREPAREGENOME.out.orfs_genomic_bed  // channel: [ val(meta), path(*.orfs-genomic.annotated.bed.gz) ]
    orfs_exons_bed   = RPBP_PREPAREGENOME.out.orfs_exons_bed    // channel: [ val(meta), path(*.orfs-exons.annotated.bed.gz) ]
    predicted        = BAM_RPBP_PREDICTORFS.out.predicted       // channel: [ val(meta), path(*.predicted-orfs.filtered.bed.gz) ]
    dna_fasta        = BAM_RPBP_PREDICTORFS.out.dna_fasta       // channel: [ val(meta), path(*.predicted-orfs.filtered.dna.fa) ]
    protein_fasta    = BAM_RPBP_PREDICTORFS.out.protein_fasta   // channel: [ val(meta), path(*.predicted-orfs.filtered.protein.fa) ]
    orf_bayes        = BAM_RPBP_PREDICTORFS.out.orf_bayes       // channel: [ val(meta), path(*.bayes-factors.bed.gz) ]
}
