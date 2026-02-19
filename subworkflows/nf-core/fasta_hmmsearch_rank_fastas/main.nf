include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMER_HMMRANK   } from '../../../modules/nf-core/hmmer/hmmrank/main'
include { SEQTK_SUBSEQ    } from '../../../modules/nf-core/seqtk/subseq/main'

workflow FASTA_HMMSEARCH_RANK_FASTAS {

    take:
    ch_hmms  // channel: [ val(meta), file(hmm) ], i.e. a list of hmm profiles, each with its meta object
    ch_fasta // channel: file(fasta), a single fasta file

    main:

    ch_versions = channel.empty()

    ch_hmms
        .combine(ch_fasta)
        .map { index -> [ index[0], index[1], index[2], false, true, false ] }
        .set { ch_hmmsearch }

    HMMER_HMMSEARCH ( ch_hmmsearch )
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions.first())

    HMMER_HMMSEARCH.out.target_summary
        .collect { index -> index[1] }
        .map { index -> [ [ id: 'rank' ], index ] }
        .set { ch_hmmrank }

    HMMER_HMMRANK ( ch_hmmrank )
    ch_versions = ch_versions.mix(HMMER_HMMRANK.out.versions.first())

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
    hmmrank                 = HMMER_HMMRANK.out.hmmrank       // channel: [ [ id: 'rank' ], hmmrank_tsv ]
    seqfastas               = SEQTK_SUBSEQ.out.sequences      // channel: [ meta, fasta ]

    versions                = ch_versions                     // channel: [ versions.yml ]
}
