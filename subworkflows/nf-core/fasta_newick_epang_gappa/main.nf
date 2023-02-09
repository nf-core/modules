include { HMMER_HMMBUILD                            } from '../../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMALIGN as HMMER_HMMALIGNREF       } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_HMMALIGN as HMMER_HMMALIGNQUERY     } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK as HMMER_MASKREF         } from '../../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLALIMASK as HMMER_MASKQUERY       } from '../../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT as HMMER_UNALIGNREF     } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as HMMER_AFAFORMATREF   } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as HMMER_AFAFORMATQUERY } from '../../../modules/nf-core/hmmer/eslreformat/main'
include { MAFFT                                     } from '../../../modules/nf-core/mafft/main'
include { EPANG_PLACE                               } from '../../../modules/nf-core/epang/place/main'
include { EPANG_SPLIT                               } from '../../../modules/nf-core/epang/split/main'
include { GAPPA_EXAMINEGRAFT as GAPPA_GRAFT         } from '../../../modules/nf-core/gappa/examinegraft/main'
include { GAPPA_EXAMINEASSIGN as GAPPA_ASSIGN       } from '../../../modules/nf-core/gappa/examineassign/main'
include { GAPPA_EXAMINEHEATTREE as GAPPA_HEATTREE   } from '../../../modules/nf-core/gappa/examineheattree/main'

workflow FASTA_NEWICK_EPANG_GAPPA {

    take:
    ch_pp_data // channel: [ meta: val(meta), data: [ alignmethod: val(alignmethod), queryseqfile: file(queryseqfile), refseqfile: file(refseqfile), refphylogeny: file(refphylogeny), hmmfile: file(hmmfile), model: val(model) ] ]

    main:
    ch_versions = Channel.empty()

    // Divide the input channel into two: One for hmmer and one for mafft alignment
    ch_hmmer_data = ch_pp_data.filter { it.data.alignmethod == 'hmmer' }
    ch_mafft_data = ch_pp_data.filter { it.data.alignmethod == 'mafft' }

    // 1.a.1 HMMER alignment: For entries that do not specify an hmm file, build one to use for alignment
    HMMER_HMMBUILD (
        ch_hmmer_data
            .filter { ! it.data.hmmfile }
            .map { [ it.meta, it.data.refseqfile ] },
        []
    )
    // 1.a.2 This handles mixed input where some samples have hmmfile set, while others don't (sample sheet input)
    ch_hmm = Channel.empty()
        .mix(HMMER_HMMBUILD.out.hmm.map { [ it[0], it[1] ] })
        .mix(
            ch_hmmer_data
                .filter { it.data.hmmfile }
                .map { [ it.meta, it.data.hmmfile ] }
        )

    ch_versions = ch_versions.mix(HMMER_HMMBUILD.out.versions.first())

    // 1.b For entries that do not specify an hmm file, "unalign" the reference sequences before they can be aligned to the hmm.
    HMMER_UNALIGNREF (
        ch_hmmer_data
            .filter { ! it.data.hmmfile }
            .map { [ it.meta, it.data.refseqfile ] }
    )
    ch_hmmer_unaligned = Channel.empty()
        .mix(HMMER_UNALIGNREF.out.seqreformated.map { [ it[0], it[1] ] })
        .mix(
            ch_hmmer_data
                .filter { it.data.hmmfile }
                .map { [ it.meta, it.data.refseqfile ] }
        )

    ch_versions = ch_versions.mix(HMMER_UNALIGNREF.out.versions)

    // 1.c Align the reference and query sequences to the profile
    ch_hmmer_alignref = ch_hmm
        .mix(ch_hmmer_unaligned)
        .groupTuple(size: 2, sort: { a, b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNREF (
        ch_hmmer_alignref.map { [ it[0], it[1][0] ] },
        ch_hmmer_alignref.map { it[1][1] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMALIGNREF.out.versions)

    ch_hmmer_alignquery = Channel.empty()
        .mix(ch_hmmer_data.map { [ it.meta, it.data.queryseqfile ] })
        .mix(ch_hmm)
        .groupTuple(size: 2, sort: { a, b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNQUERY (
        ch_hmmer_alignquery.map { [ it[0], it[1][0] ] },
        ch_hmmer_alignquery.map { it[1][1] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMALIGNQUERY.out.versions)

    // 1.d Mask the alignments (Add '--rf-is-mask' ext.args in config for the process.)
    HMMER_MASKREF ( HMMER_HMMALIGNREF.out.sthlm.map { [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_MASKREF.out.versions)

    HMMER_MASKQUERY ( HMMER_HMMALIGNQUERY.out.sthlm.map { [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_MASKQUERY.out.versions)

    // 1.e Reformat alignments to "afa" (aligned fasta)
    HMMER_AFAFORMATREF ( HMMER_MASKREF.out.maskedaln )
    ch_versions = ch_versions.mix(HMMER_AFAFORMATREF.out.versions)

    HMMER_AFAFORMATQUERY ( HMMER_MASKQUERY.out.maskedaln )
    ch_versions = ch_versions.mix(HMMER_AFAFORMATQUERY.out.versions)

    // 2.a MAFFT profile alignment of query sequences to reference alignment
    MAFFT (
        ch_mafft_data.map { [ it.meta, it.data.refseqfile ] },
        ch_mafft_data.map { [ it.data.queryseqfile ] }
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions)

    // 2.b Split the profile alignment into reference and query parts
    EPANG_SPLIT (
        ch_mafft_data.map { [ it.meta, it.data.refseqfile ] }
            .join(MAFFT.out.fas)
    )
    ch_versions = ch_versions.mix(EPANG_SPLIT.out.versions)

    // 3. Do the placement
    ch_epang_query = ch_pp_data.map { [ it.meta, it.data.model, it.data.refphylogeny ] }
        .join ( HMMER_AFAFORMATQUERY.out.seqreformated )
        .join ( HMMER_AFAFORMATREF.out.seqreformated )
        .mix(
            ch_pp_data.map { [ it.meta, it.data.model, it.data.refphylogeny ] }
                .join(EPANG_SPLIT.out.query.map { [ it[0], it[1] ] } )
                .join(EPANG_SPLIT.out.reference.map { [ it[0], it[1] ] } )
        )
        .map { [ [ id:it[0].id, model:it[1] ], it[3], it[4], it[2] ] }

    EPANG_PLACE (
        ch_epang_query,
        [], []
    )
    ch_versions = ch_versions.mix(EPANG_PLACE.out.versions)

    // 7. Calculate a tree with the placed sequences
    GAPPA_GRAFT ( EPANG_PLACE.out.jplace )
    ch_versions = ch_versions.mix(GAPPA_GRAFT.out.versions)

    // 8. Classify
    GAPPA_ASSIGN (
        EPANG_PLACE.out.jplace
            .map { [ [ id:it[0].id ], it[1] ] }
            .join( ch_pp_data.map { [ it.meta, it.data.taxonomy ] } )
    )
    ch_versions = ch_versions.mix(GAPPA_ASSIGN.out.versions)

    // 9. Heat tree output
    GAPPA_HEATTREE ( EPANG_PLACE.out.jplace )
    ch_versions = ch_versions.mix(GAPPA_HEATTREE.out.versions)

    emit:
    epang               = EPANG_PLACE.out.epang
    jplace              = EPANG_PLACE.out.jplace
    grafted_phylogeny   = GAPPA_GRAFT.out.newick
    taxonomy_profile    = GAPPA_ASSIGN.out.profile
    taxonomy_per_query  = GAPPA_ASSIGN.out.per_query
    heattree            = GAPPA_HEATTREE.out.svg
    versions            = ch_versions                     // channel: [ versions.yml ]
}

