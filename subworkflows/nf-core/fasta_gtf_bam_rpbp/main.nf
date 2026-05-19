//
// End-to-end Rp-Bp execution: render a YAML config, build the Rp-Bp genome
// index once, then run the per-sample ORF prediction chain on each Ribo-seq
// BAM. Deliberately avoids rpbp's `run-rpbp-pipeline` / `create-orf-profiles`
// wrappers - both invoke flexbar + bowtie + STAR on raw FASTQs, duplicating
// the upstream alignment supplied as input. Splitting also gives independent
// caching per step on resume.
//

include { RPBP_BUILDCONFIG       } from '../../../modules/nf-core/rpbp/buildconfig/main'
include { RPBP_PREPAREGENOME     } from '../../../modules/nf-core/rpbp/preparegenome/main'
include { BAM_RPBP_PREDICTORFS   } from '../bam_rpbp_predictorfs/main'

workflow FASTA_GTF_BAM_RPBP {

    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ] - Ribo-seq BAMs
    ch_fasta_gtf  // channel (single value): [ val(meta), path(fasta), path(gtf) ]
    extra_yaml    // val: optional extra YAML appended to the rpbp config (pass '' for defaults)

    main:

    // 1. Build the Rp-Bp YAML config from pipeline inputs.
    RPBP_BUILDCONFIG(
        ch_fasta_gtf,
        'reference',
        extra_yaml
    )

    // 2. Prepare the Rp-Bp index once per pipeline invocation.
    ch_preparegenome_in = ch_fasta_gtf
        .combine(RPBP_BUILDCONFIG.out.config.map { _meta, cfg -> cfg })
        .map { meta, fasta, gtf, cfg -> [ meta, fasta, gtf, cfg ] }

    RPBP_PREPAREGENOME(ch_preparegenome_in)

    // 3. Take dedicated emissions so downstream sees them as proper files
    //    rather than reaching into the index dir (which loses the s3:// scheme
    //    under Fusion).
    ch_transcript_bed   = RPBP_PREPAREGENOME.out.transcript_bed  .map { _meta, bed -> bed }.first()
    ch_orfs_genomic_bed = RPBP_PREPAREGENOME.out.orfs_genomic_bed.map { _meta, bed -> bed }.first()
    ch_orfs_exons_bed   = RPBP_PREPAREGENOME.out.orfs_exons_bed  .map { _meta, bed -> bed }.first()
    ch_genome_fasta     = ch_fasta_gtf.map { _meta, fasta, _gtf -> fasta }.first()

    // 4. Per-sample ORF prediction chain.
    BAM_RPBP_PREDICTORFS (
        ch_bam,
        ch_transcript_bed,
        ch_orfs_genomic_bed,
        ch_orfs_exons_bed,
        ch_genome_fasta
    )

    emit:
    config           = RPBP_BUILDCONFIG.out.config              // channel: [ val(meta), path(rpbp_config.yaml) ]
    index            = RPBP_PREPAREGENOME.out.index             // channel: [ val(meta), path(rpbp_index), path(config) ]
    transcript_bed   = RPBP_PREPAREGENOME.out.transcript_bed    // channel: [ val(meta), path(*.annotated.bed.gz) ]
    orfs_genomic_bed = RPBP_PREPAREGENOME.out.orfs_genomic_bed  // channel: [ val(meta), path(*.orfs-genomic.bed.gz) ]
    orfs_exons_bed   = RPBP_PREPAREGENOME.out.orfs_exons_bed    // channel: [ val(meta), path(*.orfs-exons.bed.gz) ]
    predicted        = BAM_RPBP_PREDICTORFS.out.predicted       // channel: [ val(meta), path(*.predicted-orfs.filtered.bed.gz) ]
    dna_fasta        = BAM_RPBP_PREDICTORFS.out.dna_fasta       // channel: [ val(meta), path(*.predicted-orfs.filtered.dna.fa) ]
    protein_fasta    = BAM_RPBP_PREDICTORFS.out.protein_fasta   // channel: [ val(meta), path(*.predicted-orfs.filtered.protein.fa) ]
    orf_bayes        = BAM_RPBP_PREDICTORFS.out.orf_bayes       // channel: [ val(meta), path(*.bayes-factors.bed.gz) ]
}
