include { HMMER_HMMBUILD as HMMER_HMMBUILD          } from '../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMALIGN as HMMER_HMMALIGNREF       } from '../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_HMMALIGN as HMMER_HMMALIGNQUERY     } from '../../modules/nf-core/hmmer/hmmalign/main'
include { HMMER_ESLALIMASK as HMMER_MASKREF         } from '../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLALIMASK as HMMER_MASKQUERY       } from '../../modules/nf-core/hmmer/eslalimask/main'
include { HMMER_ESLREFORMAT as HMMER_UNALIGNREF     } from '../../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as HMMER_AFAFORMATREF   } from '../../modules/nf-core/hmmer/eslreformat/main'
include { HMMER_ESLREFORMAT as HMMER_AFAFORMATQUERY } from '../../modules/nf-core/hmmer/eslreformat/main'
include { EPANG                                     } from '../../modules/nf-core/epang/main'
include { GAPPA_EXAMINEGRAFT as GAPPA_GRAFT         } from '../../modules/nf-core/gappa/examinegraft/main'
include { GAPPA_EXAMINEASSIGN as GAPPA_ASSIGN       } from '../../modules/nf-core/gappa/examineassign/main'
include { GAPPA_EXAMINEHEATTREE as GAPPA_HEATTREE   } from '../../modules/nf-core/gappa/examineheattree/main'

workflow SEQUENCES_NEWICK_PHYLOPLACE_EPANG_GAPPA {

    take:
    ch_pp_data // channel: [ val(meta), [ file(queryseqfile), file(refseqfile), file(refphylogeny), file(hmmfile), val(model) ] ]

    main:

    ch_versions = Channel.empty()

    // 1.a For entries that do not specify an hmm file, build one to use for alignment
    HMMER_HMMBUILD ( 
        ch_pp_data
            .filter { ! it.data.hmmfile }
            .map { [ it.meta, it.data.refseqfile ] }, 
        [] 
    )
    // 1.b This handles mixed input where some samples have hmmfile set, while others don't (sample sheet input)
    ch_hmm = Channel.empty()
        .mix(HMMER_HMMBUILD.out.hmm.map { [ it[0], it[1] ] })
        .mix(
            ch_pp_data
                .filter { it.data.hmmfile }
                .map { [ it.meta, it.data.hmmfile ] }
        )

    ch_versions = ch_versions.mix(HMMER_HMMBUILD.out.versions.first())

    // 2. For entries that do not specify an hmm file, "unalign" the reference sequences before they can be aligned to the hmm.
    HMMER_UNALIGNREF ( 
        ch_pp_data
            .filter { ! it.data.hmmfile }
            .map { [ it.meta, it.data.refseqfile ] } 
    )
    ch_unaligned = Channel.empty()
        .mix(HMMER_UNALIGNREF.out.seqreformated.map { [ it[0], it[1] ] })
        .mix(
            ch_pp_data
                .filter { it.data.hmmfile }
                .map { [ it.meta, it.data.refseqfile ] }
        )

    ch_versions = ch_versions.mix(HMMER_UNALIGNREF.out.versions)

    // 3. Align the reference and query sequences to the profile
    ch_alignref = ch_hmm
        .mix(ch_unaligned)
        .groupTuple(size: 2, sort: { a, b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNREF ( 
        ch_alignref.map { [ it[0], it[1][0] ] },
        ch_alignref.map { it[1][1] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMALIGNREF.out.versions)

    ch_alignquery = Channel.empty()
        .mix(ch_pp_data.map { [ it.meta, it.data.queryseqfile ] })
        .mix(ch_hmm)
        .groupTuple(size: 2, sort: { a, b -> a =~ /\.hmm/ ? 1 : -1 })

    HMMER_HMMALIGNQUERY (
        ch_alignquery.map { [ it[0], it[1][0] ] },
        ch_alignquery.map { it[1][1] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMALIGNQUERY.out.versions)

    // 4. Mask the alignments (Add '--rf-is-mask' ext.args in config for the process.)
    HMMER_MASKREF ( HMMER_HMMALIGNREF.out.sthlm.map { [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_MASKREF.out.versions)

    HMMER_MASKQUERY ( HMMER_HMMALIGNQUERY.out.sthlm.map { [ it[0], it[1], [], [], [], [], [], [] ] }, [] )
    ch_versions = ch_versions.mix(HMMER_MASKQUERY.out.versions)

    // 5. Reformat alignments to "afa" (aligned fasta)
    HMMER_AFAFORMATREF ( HMMER_MASKREF.out.maskedaln )
    ch_versions = ch_versions.mix(HMMER_AFAFORMATREF.out.versions)

    HMMER_AFAFORMATQUERY ( HMMER_MASKQUERY.out.maskedaln )
    ch_versions = ch_versions.mix(HMMER_AFAFORMATQUERY.out.versions)

    // 6. Do the placement
    ch_epang_query = ch_pp_data.map { [ it.meta, it.data.model ] }
        .join ( HMMER_AFAFORMATQUERY.out.seqreformated )
        .map { [ [ id:it[0].id, model:it[1] ], it[2] ] }
    EPANG (
        ch_epang_query,
        HMMER_AFAFORMATREF.out.seqreformated.map { it[1] },
        ch_pp_data.map { it.data.refphylogeny },
        [], [], []
    )
    ch_versions = ch_versions.mix(EPANG.out.versions)

    // 7. Calculate a tree with the placed sequences
    GAPPA_GRAFT ( EPANG.out.jplace )
    ch_versions = ch_versions.mix(GAPPA_GRAFT.out.versions)

    // 8. Classify
    GAPPA_ASSIGN ( EPANG.out.jplace, ch_pp_data.map { it.data.taxonomy } )
    ch_versions = ch_versions.mix(GAPPA_ASSIGN.out.versions)

    // 9. Heat tree output
    GAPPA_HEATTREE ( EPANG.out.jplace )
    ch_versions = ch_versions.mix(GAPPA_HEATTREE.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

