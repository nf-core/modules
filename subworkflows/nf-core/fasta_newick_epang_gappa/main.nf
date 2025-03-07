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
include { EPANG_PLACE                               } from '../../../modules/nf-core/epang/place/main'
include { EPANG_SPLIT as EPANG_SPLIT_CLUSTALO       } from '../../../modules/nf-core/epang/split/main'
include { EPANG_SPLIT as EPANG_SPLIT_MAFFT          } from '../../../modules/nf-core/epang/split/main'
include { GAPPA_EXAMINEGRAFT as GAPPA_GRAFT         } from '../../../modules/nf-core/gappa/examinegraft/main'
include { GAPPA_EXAMINEASSIGN as GAPPA_ASSIGN       } from '../../../modules/nf-core/gappa/examineassign/main'
include { GAPPA_EXAMINEHEATTREE as GAPPA_HEATTREE   } from '../../../modules/nf-core/gappa/examineheattree/main'

workflow FASTA_NEWICK_EPANG_GAPPA {

    take:
    ch_pp_data // channel: [ meta: val(meta), data: [ alignmethod: val(alignmethod), queryseqfile: file(queryseqfile), refseqfile: file(refseqfile), refphylogeny: file(refphylogeny), hmmfile: file(hmmfile), model: val(model) ] ]

    main:
    ch_versions = Channel.empty()

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
    ch_hmm = Channel.empty()
        .mix(HMMER_HMMBUILD.out.hmm.map { it -> [ it[0], it[1] ] })
        .mix(
            ch_hmmer_data
                .filter { it -> it.data.hmmfile }
                .map { it -> [ it.meta, it.data.hmmfile ] }
        )

    ch_versions = ch_versions.mix(HMMER_HMMBUILD.out.versions.first())

    // 1.b For entries that do not specify an hmm file, "unalign" the reference sequences before they can be aligned to the hmm.
    HMMER_UNALIGNREF (
        ch_hmmer_data
            .filter { it -> ! it.data.hmmfile }
            .map { it -> [ it.meta, it.data.refseqfile ] }
    )
    ch_hmmer_unaligned = Channel.empty()
        .mix(HMMER_UNALIGNREF.out.seqreformated.map { it -> [ it[0], it[1] ] })
        .mix(
            ch_hmmer_data
                .filter { it -> it.data.hmmfile }
                .map { it -> [ it.meta, it.data.refseqfile ] }
        )

    ch_versions = ch_versions.mix(HMMER_UNALIGNREF.out.versions)

    // 1.c Align the reference and query sequences to the profile
    ch_hmmer_alignref = ch_hmm
        .mix(ch_hmmer_unaligned)
        .groupTuple(size: 2, sort: { a, b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNREF (
        ch_hmmer_alignref.map { it -> [ it[0], it[1][0] ] },
        ch_hmmer_alignref.map { it -> it[1][1] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMALIGNREF.out.versions)

    ch_hmmer_alignquery = Channel.empty()
        .mix(ch_hmmer_data.map { it -> [ it.meta, it.data.queryseqfile ] })
        .mix(ch_hmm)
        .groupTuple(size: 2, sort: { a, b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNQUERY (
        ch_hmmer_alignquery.map { it -> [ it[0], it[1][0] ] },
        ch_hmmer_alignquery.map { it -> it[1][1] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMALIGNQUERY.out.versions)

    // 1.d Mask the alignments (Add '--rf-is-mask' ext.args in config for the process.)
    HMMER_MASKREF ( HMMER_HMMALIGNREF.out.sto.map { it -> [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_MASKREF.out.versions)

    HMMER_MASKQUERY ( HMMER_HMMALIGNQUERY.out.sto.map { it -> [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_MASKQUERY.out.versions)

    // 1.e Reformat alignments to "afa" (aligned fasta)
    HMMER_AFAFORMATREF ( HMMER_MASKREF.out.maskedaln )
    ch_versions = ch_versions.mix(HMMER_AFAFORMATREF.out.versions)

    HMMER_AFAFORMATQUERY ( HMMER_MASKQUERY.out.maskedaln )
    ch_versions = ch_versions.mix(HMMER_AFAFORMATQUERY.out.versions)

    // 2.a CLUSTALO_ALIGN profile alignment of query sequences to reference alignment
    CLUSTALO_ALIGN (
        ch_clustalo_data.map { it -> [ it.meta, it.data.queryseqfile ] },
        [ [:], []],
        [ ],
        [ ],
        ch_clustalo_data.map { it -> it.data.refseqfile },
        [ ],
        false
    )
    ch_versions = ch_versions.mix(CLUSTALO_ALIGN.out.versions)

    // 2.b Split the profile alignment into reference and query parts
    EPANG_SPLIT_CLUSTALO (
        ch_clustalo_data.map { it -> [ it.meta, it.data.refseqfile ] }
            .join(CLUSTALO_ALIGN.out.alignment)
    )
    ch_versions = ch_versions.mix(EPANG_SPLIT_CLUSTALO.out.versions)

    // 3.a MAFFT profile alignment of query sequences to reference alignment
    MAFFT_ALIGN (
        ch_mafft_data.map { it -> [ it.meta, it.data.refseqfile ] },
        ch_mafft_data.map { it -> [ it.meta, it.data.queryseqfile ] },
        [ [], [] ],
        [ [], [] ],
        [ [], [] ],
        [ [], [] ],
        false
    )
    ch_versions = ch_versions.mix(MAFFT_ALIGN.out.versions)

    // 3.b Split the profile alignment into reference and query parts
    EPANG_SPLIT_MAFFT (
        ch_mafft_data.map { it -> [ it.meta, it.data.refseqfile ] }
            .join(MAFFT_ALIGN.out.fas)
    )
    ch_versions = ch_versions.mix(EPANG_SPLIT_MAFFT.out.versions)

    // 4. Do the placement
    ch_epang_query = ch_pp_data.map { it -> [ it.meta, it.data.model, it.data.refphylogeny ] }
        .join ( HMMER_AFAFORMATQUERY.out.seqreformated )
        .join ( HMMER_AFAFORMATREF.out.seqreformated )
        .mix(
            ch_pp_data.map { it -> [ it.meta, it.data.model, it.data.refphylogeny ] }
                .join(
                    EPANG_SPLIT_CLUSTALO.out.query
                        .mix(EPANG_SPLIT_MAFFT.out.query)
                        .map { it -> [ it[0], it[1] ] }
                )
                .join(
                    EPANG_SPLIT_CLUSTALO.out.reference
                        .mix(EPANG_SPLIT_MAFFT.out.reference)
                        .map { it -> [ it[0], it[1] ] }
                )
        )
        .map { it -> [ [ id:it[0].id, model:it[1] ], it[3], it[4], it[2] ] }

    EPANG_PLACE (
        ch_epang_query,
        [], []
    )
    ch_versions = ch_versions.mix(EPANG_PLACE.out.versions)

    // 5. Calculate a tree with the placed sequences
    GAPPA_GRAFT ( EPANG_PLACE.out.jplace )
    ch_versions = ch_versions.mix(GAPPA_GRAFT.out.versions)

    // 6. Classify
    GAPPA_ASSIGN (
        EPANG_PLACE.out.jplace
            .map { it -> [ [ id:it[0].id ], it[1] ] }
            .join( ch_pp_data.map { it -> [ [ id: it.meta.id ], it.data.taxonomy ] } )
    )
    ch_versions = ch_versions.mix(GAPPA_ASSIGN.out.versions)

    // 7. Heat tree output
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
