include { CUSTOM_SHORTENFASTAIDS    } from '../../../modules/nf-core/custom/shortenfastaids/main'
include { LTRHARVEST                } from '../../../modules/nf-core/ltrharvest/main'
include { LTRFINDER                 } from '../../../modules/nf-core/ltrfinder/main'
include { LTRRETRIEVER_LTRRETRIEVER } from '../../../modules/nf-core/ltrretriever/ltrretriever/main'
include { CAT_CAT                   } from '../../../modules/nf-core/cat/cat/main'
include { LTRRETRIEVER_LAI          } from '../../../modules/nf-core/ltrretriever/lai/main'
include { CUSTOM_RESTOREGFFIDS      } from '../../../modules/nf-core/custom/restoregffids/main'

workflow FASTA_LTRRETRIEVER_LAI {

    take:
    ch_fasta                        // channel: [ val(meta), fasta ]
    ch_monoploid_seqs               // channel: [ val(meta), txt ]; Optional: Set to [] if not needed
                                    // val(meta) from ch_fasta and ch_monoploid_seqs are only required
                                    // to have the same `id`
    skip_lai                        // val(true|false)

    main:
    ch_versions                     = Channel.empty()

    // MOUDLE: CUSTOM_SHORTENFASTAIDS
    CUSTOM_SHORTENFASTAIDS ( ch_fasta )

    ch_short_ids_fasta              = ch_fasta
                                    | join(CUSTOM_SHORTENFASTAIDS.out.short_ids_fasta, by:0, remainder:true)
                                    | map { meta, fasta, short_ids_fasta ->
                                        if ( fasta ) { [ meta, short_ids_fasta ?: fasta ] }
                                    }

    ch_short_ids_tsv                = CUSTOM_SHORTENFASTAIDS.out.short_ids_tsv
    ch_short_monoploid_seqs         = ch_short_ids_tsv
                                    | join(
                                        ch_monoploid_seqs ?: Channel.empty()
                                    )
                                    | filter { meta, tsv, seqs -> seqs } // Cater to channel: [ meta, tsv, [] ]
                                    | map { meta, short_ids_tsv, monoploid_seqs ->
                                        map_monoploid_seqs_to_new_ids(meta, short_ids_tsv, monoploid_seqs)
                                    }
                                    | collectFile(newLine:true)
                                    | map { seqs ->
                                        def id = seqs.name.split('.mapped.monoploid.seqs.txt')[0]

                                        [ [ id: id ], seqs ]
                                    }
    ch_versions                     = ch_versions.mix(CUSTOM_SHORTENFASTAIDS.out.versions.first())

    // MODULE: LTRHARVEST
    LTRHARVEST ( ch_short_ids_fasta )

    ch_ltrharvest_scn               = LTRHARVEST.out.scn
    ch_versions                     = ch_versions.mix(LTRHARVEST.out.versions.first())

    // MODULE: LTRFINDER
    LTRFINDER ( ch_short_ids_fasta )

    ch_ltrfinder_scn                = LTRFINDER.out.scn
    ch_versions                     = ch_versions.mix(LTRFINDER.out.versions.first())

    // MODULE: CAT_CAT
    ch_cat_cat_inputs               = ch_ltrharvest_scn
                                    | join(ch_ltrfinder_scn)
                                    | map { meta, harvested, found -> [ meta, [ harvested, found ] ] }

    CAT_CAT ( ch_cat_cat_inputs )

    ch_ltr_candidates               = CAT_CAT.out.file_out
    ch_versions                     = ch_versions.mix(CAT_CAT.out.versions.first())

    // MODULE: LTRRETRIEVER_LTRRETRIEVER
    ch_ltrretriever_inputs          = ch_short_ids_fasta.join(ch_ltr_candidates)

    LTRRETRIEVER_LTRRETRIEVER (
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_pass_list                    = LTRRETRIEVER_LTRRETRIEVER.out.pass_list
    ch_ltrlib                       = LTRRETRIEVER_LTRRETRIEVER.out.ltrlib
    ch_annotation_out               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_out
    ch_annotation_gff               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: LAI
    ch_lai_inputs                   = skip_lai
                                    ? Channel.empty()
                                    : ch_short_ids_fasta
                                    | join(ch_pass_list)
                                    | join(ch_annotation_out)
                                    | map { meta, fasta, pass, out ->
                                        [ meta.id, meta, fasta, pass, out ]
                                    }
                                    | join(
                                        ch_short_monoploid_seqs
                                        | map { meta, mono -> [ meta.id, mono ] },
                                        by:0,
                                        remainder: true
                                    )
                                    | map { id, meta, fasta, pass, out, mono ->
                                        [ meta, fasta, pass, out, mono ?: [] ]
                                    }
    LTRRETRIEVER_LAI(
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> [ meta, fasta ] },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> pass },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> out },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> mono }
    )

    ch_lai_log                      = LTRRETRIEVER_LAI.out.log
    ch_lai_out                      = LTRRETRIEVER_LAI.out.lai_out
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LAI.out.versions.first())

    // MODULE: CUSTOM_RESTOREGFFIDS
    ch_restorable_gff_tsv           = ch_annotation_gff.join(ch_short_ids_tsv)

    CUSTOM_RESTOREGFFIDS (
        ch_restorable_gff_tsv.map { meta, gff, tsv -> [ meta, gff ] },
        ch_restorable_gff_tsv.map { meta, gff, tsv -> tsv }
    )

    ch_restored_gff                 = ch_annotation_gff
                                    | join(CUSTOM_RESTOREGFFIDS.out.restored_ids_gff3, by:0, remainder:true)
                                    | map { meta, gff, restored_gff -> [ meta, restored_gff ?: gff ] }

    ch_versions                     = ch_versions.mix(CUSTOM_RESTOREGFFIDS.out.versions.first())

    emit:
    ltrlib                          = ch_ltrlib         // channel: [ val(meta), fasta ]
    annotation_gff                  = ch_restored_gff   // channel: [ val(meta), gff ]
    lai_log                         = ch_lai_log        // channel: [ val(meta), log ]
    lai_out                         = ch_lai_out        // channel: [ val(meta), out ]
    versions                        = ch_versions       // channel: [ versions.yml ]
}


def map_monoploid_seqs_to_new_ids(meta, short_ids_tsv, monoploid_seqs) {

    def short_ids_head = short_ids_tsv.text.split('\n')[0]

    if (short_ids_head == "IDs have acceptable length and character. No change required.") {
        return [ "${meta.id}.mapped.monoploid.seqs.txt" ] + monoploid_seqs.text.split('\n')
    }

    def orig_to_new_ids = [:]
    short_ids_tsv.text.eachLine { line ->
        def (original_id, renamed_id) = line.split('\t')
        orig_to_new_ids[original_id] = renamed_id
    }

    def mapped_ids = []
    monoploid_seqs.text.eachLine { original_id ->
        if (!orig_to_new_ids[original_id]) {
            error "Faild to find $original_id in ${monoploid_seqs}" +
            "The monoploid_seqs file is malformed!"
        }

        mapped_ids.add(orig_to_new_ids[original_id])
    }

    return [ "${meta.id}.mapped.monoploid.seqs.txt" ] + mapped_ids
}
