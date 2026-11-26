// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

include { HMMER_HMMBUILD                            } from '../../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMALIGN as HMMER_HMMALIGNREF       } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_HMMALIGN as HMMER_HMMALIGNQUERY     } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK as HMMER_MASKREF         } from '../../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLALIMASK as HMMER_MASKQUERY       } from '../../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT as HMMER_UNALIGNREF     } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as HMMER_AFAFORMATREF   } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as HMMER_AFAFORMATQUERY } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { CLUSTALO_ALIGN                            } from '../../../modules/nf-core/clustalo/align/main'
include { MAFFT_ALIGN                               } from '../../../modules/nf-core/mafft/align/main'
include { EPANG_FILTERQUERY                         } from '../../../modules/nf-core/epang/filterquery/main'
include { EPANG_PLACE                               } from '../../../modules/nf-core/epang/place/main'
include { EPANG_SPLIT as EPANG_SPLIT_CLUSTALO       } from '../../../modules/nf-core/epang/split/main'
include { EPANG_SPLIT as EPANG_SPLIT_MAFFT          } from '../../../modules/nf-core/epang/split/main'
include { GAPPA_EXAMINEGRAFT as GAPPA_GRAFT         } from '../../../modules/nf-core/gappa/examinegraft/main'
include { GAPPA_EXAMINEASSIGN as GAPPA_ASSIGN       } from '../../../modules/nf-core/gappa/examineassign/main'
include { GAPPA_EXAMINEHEATTREE as GAPPA_HEATTREE   } from '../../../modules/nf-core/gappa/examineheattree/main'

workflow FASTA_NEWICK_EPANG_GAPPA {

    take:
    ch_pp_data // channel: [ meta: val(meta), data: [ alignmethod: val(alignmethod), queryseqfile: file(queryseqfile), refseqfile: file(refseqfile), refphylogeny: file(refphylogeny), hmmfile: file(hmmfile), model: val(model) ] ]
    compress_alignment // value: boolean, write the clustalo and mafft alignments gzipped

    main:

    // Divide the input channel into three: One each for hmmer, clustalo and mafft alignment
    ch_hmmer_data    = ch_pp_data.filter { it -> it.data.alignmethod == 'hmmer'    }
    ch_clustalo_data = ch_pp_data.filter { it -> it.data.alignmethod == 'clustalo' }
    ch_mafft_data    = ch_pp_data.filter { it -> it.data.alignmethod == 'mafft'    }

    // 1.a.1 HMMER alignment: For entries that do not specify an hmm file, build one to use for alignment
    HMMER_HMMBUILD (
        ch_hmmer_data
            .filter { it -> ! it.data.hmmfile }
            .map { it -> [ it.meta, it.data.refseqfile ] },
        []
    )
    // 1.a.2 This handles mixed input where some samples have hmmfile set, while others don't (sample sheet input)
    ch_hmm = channel.empty()
        .mix(HMMER_HMMBUILD.out.hmm.map { it -> [ it[0], it[1] ] })
        .mix(
            ch_hmmer_data
                .filter { it -> it.data.hmmfile }
                .map { it -> [ it.meta, it.data.hmmfile ] }
        )

    // 1.b For entries that do not specify an hmm file, "unalign" the reference sequences before they can be aligned to the hmm.
    HMMER_UNALIGNREF (
        ch_hmmer_data
            .filter { it -> ! it.data.hmmfile }
            .map { it -> [ it.meta, it.data.refseqfile ] },
        '| sed "/^>/!s/-//g"'
    )
    ch_hmmer_unaligned = channel.empty()
        .mix(HMMER_UNALIGNREF.out.seqreformated.map { it -> [ it[0], it[1] ] })
        .mix(
            ch_hmmer_data
                .filter { it -> it.data.hmmfile }
                .map { it -> [ it.meta, it.data.refseqfile ] }
        )

    // 1.c Align the reference and query sequences to the profile
    ch_hmmer_alignref = ch_hmm
        .mix(ch_hmmer_unaligned)
        .groupTuple(size: 2, sort: { a, _b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNREF (
        ch_hmmer_alignref.map { it -> [ it[0], it[1][0] ] },
        ch_hmmer_alignref.map { it -> it[1][1] }
    )

    ch_hmmer_alignquery = channel.empty()
        .mix(ch_hmmer_data.map { it -> [ it.meta, it.data.queryseqfile ] })
        .mix(ch_hmm)
        .groupTuple(size: 2, sort: { a, _b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNQUERY (
        ch_hmmer_alignquery.map { it -> [ it[0], it[1][0] ] },
        ch_hmmer_alignquery.map { it -> it[1][1] }
    )

    // 1.d Mask the alignments (Add '--rf-is-mask' ext.args in config for the process.)
    HMMER_MASKREF ( HMMER_HMMALIGNREF.out.sto.map { it -> [ it[0], it[1], [], [], [], [], [], [] ] }, [] )

    HMMER_MASKQUERY ( HMMER_HMMALIGNQUERY.out.sto.map { it -> [ it[0], it[1], [], [], [], [], [], [] ] }, [] )

    // 1.e Reformat alignments to "afa" (aligned fasta)
    HMMER_AFAFORMATREF ( HMMER_MASKREF.out.maskedaln, '' )

    HMMER_AFAFORMATQUERY ( HMMER_MASKQUERY.out.maskedaln, '' )

    // 2.a CLUSTALO_ALIGN profile alignment of query sequences to reference alignment
    CLUSTALO_ALIGN (
        ch_clustalo_data.map { it -> [ it.meta, it.data.queryseqfile ] },
        [ [:], []],
        [ ],
        [ ],
        ch_clustalo_data.map { it -> it.data.refseqfile },
        [ ],
        compress_alignment
    )

    // 2.b Split the profile alignment into reference and query parts
    EPANG_SPLIT_CLUSTALO (
        ch_clustalo_data.map { it -> [ it.meta, it.data.refseqfile ] }
            .join(CLUSTALO_ALIGN.out.alignment)
    )

    // 3.a MAFFT profile alignment of query sequences to reference alignment
    MAFFT_ALIGN (
        ch_mafft_data.map { it -> [ it.meta, it.data.refseqfile ] },
        ch_mafft_data.map { it -> [ it.meta, it.data.queryseqfile ] },
        [ [], [] ],
        [ [], [] ],
        [ [], [] ],
        [ [], [] ],
        compress_alignment
    )

    // 3.b Split the profile alignment into reference and query parts
    EPANG_SPLIT_MAFFT (
        ch_mafft_data.map { it -> [ it.meta, it.data.refseqfile ] }
            .join(MAFFT_ALIGN.out.fas)
    )

    // 4. Prepare the placement
    ch_epang_hmmer = ch_hmmer_data
        .map { it -> [ it.meta, it.data.model, it.data.refphylogeny ] }
        .join ( HMMER_AFAFORMATQUERY.out.seqreformated )
        .join ( HMMER_AFAFORMATREF.out.seqreformated )

    ch_epang_clustalo = ch_clustalo_data
        .map { it -> [ it.meta, it.data.model, it.data.refphylogeny ] }
        .join ( EPANG_SPLIT_CLUSTALO.out.query )
        .join ( EPANG_SPLIT_CLUSTALO.out.reference )

    ch_epang_mafft = ch_pp_data
        .filter { it -> it.data.alignmethod == 'mafft' }
        .map { it -> [ it.meta, it.data.model, it.data.refphylogeny ] }
        .join ( EPANG_SPLIT_MAFFT.out.query )
        .join ( EPANG_SPLIT_MAFFT.out.reference )

    ch_epang_aligned = ch_epang_hmmer
        .mix(ch_epang_clustalo)
        .mix(ch_epang_mafft)

    // 4.a Drop the queries EPA-ng cannot place, i.e. those with no residue in the columns the
    // reference occupies: EPA-ng aborts the whole run on such a query. The thresholds are
    // optional keys of the data map, so callers that don't set them get the defaults. A sample
    // with no query left has no EPANG_FILTERQUERY query output and so no placement outputs.
    EPANG_FILTERQUERY (
        ch_epang_aligned
            .join ( ch_pp_data.map { it -> [ it.meta, it.data.min_query_columns != null ? it.data.min_query_columns : 1, it.data.min_reference_coverage != null ? it.data.min_reference_coverage : 0 ] } )
            .map { meta, _model, _tree, query, reference, min_query_columns, min_reference_coverage ->
                [ meta, reference, query, min_query_columns, min_reference_coverage ]
            }
    )

    EPANG_FILTERQUERY.out.table.subscribe { meta, table ->
        def rows    = table.readLines().drop(1)
        def dropped = rows.count { row -> row.endsWith('\tfalse') }
        if ( rows.size() > 0 && dropped == rows.size() ) {
            log.warn "EPANG_FILTERQUERY: none of the ${rows.size()} query sequences of '${meta.id}' passed the filter, so nothing is placed; see ${table.name}"
        } else if ( dropped > 0 ) {
            log.warn "EPANG_FILTERQUERY: ${dropped} of ${rows.size()} query sequences of '${meta.id}' were dropped before placement because of too few residues in the occupied reference columns; see ${table.name}"
        }
    }

    // 4.b Do the placement with the queries that passed
    ch_epang_query = ch_epang_aligned
        .map { it -> [ it[0], it[1], it[2], it[4] ] }
        .join ( EPANG_FILTERQUERY.out.query )
        .map { it -> [ [ id:it[0].id, model:it[1] ], it[4], it[3], it[2] ] }

    EPANG_PLACE (
        ch_epang_query,
        [], []
    )

    // 5. Calculate a tree with the placed sequences
    GAPPA_GRAFT ( EPANG_PLACE.out.jplace )

    // 6. Classify
    GAPPA_ASSIGN (
        EPANG_PLACE.out.jplace
            .map { it -> [ [ id:it[0].id ], it[1] ] }
            .join( ch_pp_data.map { it -> [ [ id: it.meta.id ], it.data.taxonomy ] } )
    )

    // 7. Heat tree output
    GAPPA_HEATTREE ( EPANG_PLACE.out.jplace )

    emit:
    epang               = EPANG_PLACE.out.epang
    jplace              = EPANG_PLACE.out.jplace
    epang_log           = EPANG_PLACE.out.log
    filterquery_table   = EPANG_FILTERQUERY.out.table
    grafted_phylogeny   = GAPPA_GRAFT.out.newick
    taxonomy_profile    = GAPPA_ASSIGN.out.profile
    taxonomy_per_query  = GAPPA_ASSIGN.out.per_query
    heattree            = GAPPA_HEATTREE.out.svg
    hmmbuild_log        = HMMER_HMMBUILD.out.hmmbuildout
}
