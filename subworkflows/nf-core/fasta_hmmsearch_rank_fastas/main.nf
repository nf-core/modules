include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMER_FORMATTSV as HMMER_FORMATTSV_TBLOUT       } from '../../../modules/nf-core/hmmer/formattsv/main'
include { HMMER_FORMATTSV as HMMER_FORMATTSV_DOMTBLOUT    } from '../../../modules/nf-core/hmmer/formattsv/main'
include { DUCKDB_TABLE2PARQUET as DUCKDB_TABLE2PARQUET_TBLOUT    } from '../../../modules/nf-core/duckdb/table2parquet/main'
include { DUCKDB_TABLE2PARQUET as DUCKDB_TABLE2PARQUET_DOMTBLOUT } from '../../../modules/nf-core/duckdb/table2parquet/main'
include { HMMER_HMMRANK   } from '../../../modules/nf-core/hmmer/hmmrank/main'
include { SEQTK_SUBSEQ    } from '../../../modules/nf-core/seqtk/subseq/main'

workflow FASTA_HMMSEARCH_RANK_FASTAS {

    take:
    ch_hmms         // channel: [ val(meta), file(hmm) ], i.e. a list of hmm profiles, each with its meta object
    ch_fasta        // channel: file(fasta), a single fasta file
    save_domtblout  // boolean: also write and emit hmmsearch's per-domain hit table (--domtblout)

    main:

    ch_hmms
        .combine(ch_fasta)
        .map { index -> [ index[0], index[1], index[2], false, true, save_domtblout ] }
        .set { ch_hmmsearch }

    HMMER_HMMSEARCH ( ch_hmmsearch )

    // hmmer/hmmrank ranks an already-combined, already-typed Parquet table rather than raw
    // HMMER text; hmmer/formattsv combines every profile's own table into one (labelled by
    // profile), and duckdb/table2parquet converts it. tblout and domtblout go through separate
    // calls with distinct meta ids so their Parquet filenames don't collide once both are
    // staged into the same HMMER_HMMRANK task.
    // Sorted by id: .collect() completion order isn't deterministic, and only hmmrank's own
    // output gets a final ORDER BY -- this intermediate one is published as-is.
    HMMER_HMMSEARCH.out.target_summary
        .map { meta, tbl -> [ meta.id, tbl ] }
        .collect(flat: false)
        .map { pairs -> pairs.sort { it[0] } }
        .map { pairs -> [ [ id: 'rank.tblout' ], pairs.collect { it[0] }, pairs.collect { it[1] } ] }
        .set { ch_formattsv_tblout }

    HMMER_FORMATTSV_TBLOUT ( ch_formattsv_tblout, 'tblout' )
    DUCKDB_TABLE2PARQUET_TBLOUT ( HMMER_FORMATTSV_TBLOUT.out.tsv )

    // The per-domain tables are what carry alignment coordinates, so hand them to the ranking
    // step when hmmsearch was asked to write them. save_domtblout is a plain boolean known at
    // workflow-composition time, so this if/else picks which channels get built rather than
    // branching a dataflow at runtime; ch_domtblout_parquet ends up either a single real path
    // or a single `[]`, matching what HMMER_HMMRANK treats as "no domtblout".
    if (save_domtblout) {
        // Sorted for the same reason as the tblout branch above.
        HMMER_HMMSEARCH.out.domain_summary
            .map { meta, domtbl -> [ meta.id, domtbl ] }
            .collect(flat: false)
            .map { pairs -> pairs.sort { it[0] } }
            .map { pairs -> [ [ id: 'rank.domtblout' ], pairs.collect { it[0] }, pairs.collect { it[1] } ] }
            .set { ch_formattsv_domtblout }

        HMMER_FORMATTSV_DOMTBLOUT ( ch_formattsv_domtblout, 'domtblout' )
        DUCKDB_TABLE2PARQUET_DOMTBLOUT ( HMMER_FORMATTSV_DOMTBLOUT.out.tsv )

        ch_domtblout_parquet = DUCKDB_TABLE2PARQUET_DOMTBLOUT.out.parquet.map { meta, parquet -> parquet }
    } else {
        ch_domtblout_parquet = Channel.value([])
    }

    // ch_domtblout_parquet's value is combined wrapped in an extra list, same as
    // ch_domtblouts was before this rework: combine() would otherwise treat a bare `[]`
    // value as zero items rather than one item whose payload happens to be an empty list.
    DUCKDB_TABLE2PARQUET_TBLOUT.out.parquet
        .map { meta, parquet -> [ [ id: 'rank' ], parquet ] }
        .combine(ch_domtblout_parquet.map { index -> [ index ] })
        .set { ch_hmmrank }

    HMMER_HMMRANK ( ch_hmmrank )

    HMMER_HMMRANK.out.hmmrank
        .map { index -> index[1] }
        .splitCsv(header: true, sep: '\t')
        .filter { index -> index.rank == '1' }
        .collectFile { index -> [ "${index.profile}.txt", "${index.accno}\n" ] }
        .map { index -> [ [ id: index.baseName ], index ] }
        .groupTuple(sort: true)
        .set { ch_subseq_filter }

    ch_subseq_filter
        .combine(ch_fasta)
        .map { index -> [ index[0], index[2] ] }
        .groupTuple(sort: true)
        .set { ch_subseq_fasta }

    SEQTK_SUBSEQ ( ch_subseq_fasta, ch_subseq_filter.map { index -> index[1] } )
    // SEQTK_SUBSEQ emits version as a topic channel

    emit:
    hmmrank                 = HMMER_HMMRANK.out.hmmrank             // channel: [ [ id: 'rank' ], hmmrank_tsv ]
    seqfastas               = SEQTK_SUBSEQ.out.sequences            // channel: [ meta, fasta ]
    domain_summary          = HMMER_HMMSEARCH.out.domain_summary    // channel: [ meta, domtbl ]
}
