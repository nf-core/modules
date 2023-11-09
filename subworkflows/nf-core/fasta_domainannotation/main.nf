include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTP      } from '../../../modules/nf-core/blast/blastp/main'
include { DIAMOND_MAKEDB    } from '../../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP    } from '../../../modules/nf-core/diamond/blastp/main'
include { INTERPROSCAN      } from '../../../modules/nf-core/interproscan/main'

workflow FASTA_DOMAINANNOTATION {

    take:
    ch_fasta         // channel: [ val(meta), path(fasta) ]
    val_blastp_fasta // value: /path/to/reference/fasta for blastp
    val_blastp_mode  // value: blast or diamond

    main:

    ch_versions = Channel.empty()

    if (val_blastp_mode == "blast") {
        BLAST_MAKEBLASTDB ( val_blastp_fasta )
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
        BLAST_BLASTP ( ch_fasta, BLAST_MAKEBLASTDB.out.db.first(), 'tsv' ) // .first() to also work with chunked input
        ch_versions = ch_versions.mix(BLAST_BLASTP.out.versions)
        blastp_tsv = BLAST_BLASTP.out.tsv
    } else if (val_blastp_mode == "diamond") {
        DIAMOND_MAKEDB ( val_blastp_fasta )
        ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
        blast_columns = '' // defaults: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
        DIAMOND_BLASTP ( ch_fasta, DIAMOND_MAKEDB.out.db.first(), 'txt', blast_columns )
        ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)
        blastp_tsv = DIAMOND_BLASTP.out.txt
    } else {
        throw new Exception("Invalid mode value '$val_blastp_mode'. Should be 'blast' or 'diamond'.")
    }

    INTERPROSCAN ( ch_fasta, 'tsv' )
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)

    emit:
    blastp_tsv      = blastp_tsv // channel: [ val(meta), [ tsv|txt ] ]
    inteproscan_tsv = INTERPROSCAN.out.tsv // channel: [ val(meta), [ tsv ] ]
    versions        = ch_versions // channel: [ versions.yml ]
}
