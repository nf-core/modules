include { LTRHARVEST                } from '../../../modules/nf-core/ltrharvest/main'
include { LTRFINDER                 } from '../../../modules/nf-core/ltrfinder/main'
include { LTRRETRIEVER_LTRRETRIEVER } from '../../../modules/nf-core/ltrretriever/ltrretriever/main'
include { CAT_CAT                   } from '../../../modules/nf-core/cat/cat/main'
include { LTRRETRIEVER_LAI          } from '../../../modules/nf-core/ltrretriever/lai/main'

workflow FASTA_LTRHARVEST_LTRFINDER_LTRRETRIEVER_LAI {

    take:
    ch_fasta                        // channel: [ val(meta), fasta ]
    ch_monoploid_seqs               // channel: [ val(meta2), txt ]; Optional: Set to [] if not needed
                                    // val(meta) from ch_fasta and val(meta2) from ch_monoploid_seqs are
                                    // only required to have the same `id`
    skip_lai                        // val(true|false)

    main:
    ch_versions                     = Channel.empty()

    // Prapre input channels
    ch_monoploid_seqs_plain         = ( ch_monoploid_seqs ?: Channel.empty() )
                                    | filter { meta2, seqs -> seqs }
                                    // Cater to channel: [ meta2, [] ]
                                    | map { meta2, seqs -> [ meta2.id, seqs ] }


    // collectFile: Shorten IDs
    ch_shortidstsv_orig_fasta       = ch_fasta
                                    | map { meta, fasta -> shorten_ids ( meta, fasta ) }
                                    | collectFile
                                    | map { tsv ->
                                        [ tsv.baseName.replace('.short.ids', ''), tsv ]
                                    }
                                    | join(
                                        ch_fasta.map { meta, fasta -> [ meta.id, meta, fasta ] }
                                    )
                                    | map { id, tsv, meta, fasta -> [ meta, tsv, fasta ] }

    ch_short_ids_tsv                = ch_shortidstsv_orig_fasta
                                    | map {meta, tsv, fasta -> [ meta, tsv ] }

    ch_shortidstsv_orig_fasta_branch= ch_shortidstsv_orig_fasta
                                    | branch { meta, tsv, fasta ->
                                        change: ! tsv.text.contains('IDs have acceptable length and character')
                                        nonchange: tsv.text.contains('IDs have acceptable length and character')
                                    }

    // collectFile: Create Fasta file with short IDs
    ch_short_ids_fasta              = ch_shortidstsv_orig_fasta_branch.change
                                    | map { meta, tsv, fasta ->
                                        def orig_to_new_ids = tsv.text.tokenize('\n').collectEntries { [ it.tokenize('\t')[0], it.tokenize('\t')[1] ] }

                                        def fasta_records = fasta.splitFasta( record: [ id: true, sequence: true ] )
                                        def new_fasta = fasta_records.collect { ">${orig_to_new_ids[it.id]}\n${it.sequence}" }.join('')

                                        [ "${meta.id}.short.ids.fasta", new_fasta ]
                                    }
                                    | collectFile
                                    | map { fasta -> [ fasta.baseName.replace('.short.ids', ''), fasta ] }
                                    | join(
                                        ch_shortidstsv_orig_fasta_branch.change
                                        | map { meta, tsv, fasta -> [ meta.id, meta ] }
                                    )
                                    | map { id, fasta, meta -> [ meta, fasta ] }
                                    | mix(
                                        ch_shortidstsv_orig_fasta_branch.nonchange.map { meta, tsv, fasta -> [ meta, fasta ] }
                                    )

    // collectFile: Map monoploid seqs to short IDs
    ch_short_monoploid_seqs         = ch_short_ids_tsv
                                    | map { meta, tsv -> [ meta.id, tsv ] }
                                    | join(ch_monoploid_seqs_plain)
                                    | map { id, tsv, seqs ->
                                        map_monoploid_seqs_to_new_ids(id, tsv, seqs)
                                    }
                                    | collectFile(newLine:true)
                                    | map { seqs ->
                                        def id = seqs.name.split('.mapped.monoploid.seqs.txt')[0]

                                        [ id, seqs ]
                                    }
                                    | join(
                                        ch_short_ids_tsv
                                        | map { meta, tsv -> [ meta.id, meta ] }
                                    )
                                    | map { id, seqs, meta -> [ meta, seqs ] }


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
    ch_annotation_out               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_out
    ch_pass_out                     = ch_pass_list.join(ch_annotation_out)
    ch_annotation_gff               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff
    ch_ltrlib                       = LTRRETRIEVER_LTRRETRIEVER.out.ltrlib
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: LTRRETRIEVER_LAI
    ch_short_ids_fasta_mono         = ch_short_ids_fasta
                                    | join(
                                        ch_short_monoploid_seqs,
                                        by:0,
                                        remainder: true
                                    )
                                    // Danger! This partial join can fail
                                    | filter { meta, fasta, seqs -> fasta }
                                    // This filter safeguards against fail on upstream
                                    // process failure: https://github.com/nextflow-io/nextflow/issues/5043
                                    // fasta may come from upstream processes
                                    // seqs also comes from upstream processes, it is optional
                                    // and may not be present for some of the combinations
                                    | map { meta, fasta, seqs -> [ meta, fasta, seqs ?: [] ] }

    ch_lai_inputs                   = skip_lai
                                    ? Channel.empty()
                                    : ch_short_ids_fasta_mono
                                    | join(
                                        ch_pass_out
                                    )
                                    | map { meta, fasta, seqs, pass, out ->
                                        [ meta, fasta, pass, out, seqs ]
                                    }
    LTRRETRIEVER_LAI(
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> [ meta, fasta ] },
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> pass },
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> out },
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> seqs }
    )

    ch_lai_log                      = LTRRETRIEVER_LAI.out.log
    ch_lai_out                      = LTRRETRIEVER_LAI.out.lai_out
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LAI.out.versions.first())

    // collectFile: Restore gff ids
    ch_gff_tsv_branch               = ch_annotation_gff.join(ch_short_ids_tsv)
                                    | branch { meta, gff, tsv ->
                                        change: ! tsv.text.contains('IDs have acceptable length and character')
                                        nochange: tsv.text.contains('IDs have acceptable length and character')
                                    }

    ch_collected_gff                = ch_gff_tsv_branch.change
                                    | map { meta, gff, tsv ->
                                        restore_gff_ids ( meta.id, gff, tsv )
                                    }
                                    | collectFile
                                    | map { gff ->
                                        def id = gff.name.split('.restored.ids.gff3')[0]

                                        [ id, gff ]
                                    }
                                    | join(
                                        ch_gff_tsv_branch.change
                                        | map { meta, gff, tsv -> [ meta.id, meta ] }
                                    )
                                    | map { id, gff, meta -> [ meta, gff ] }


    ch_restored_gff                 = ch_gff_tsv_branch.nochange
                                    | map { meta, gff, tsv -> [ meta, gff ] }
                                    | mix(ch_collected_gff)

    emit:
    ltrlib                          = ch_ltrlib         // channel: [ val(meta), fasta ]
    annotation_gff                  = ch_restored_gff   // channel: [ val(meta), gff ]
    lai_log                         = ch_lai_log        // channel: [ val(meta), log ]
    lai_out                         = ch_lai_out        // channel: [ val(meta), out ]
    versions                        = ch_versions       // channel: [ versions.yml ]
}


def map_monoploid_seqs_to_new_ids(id, short_ids_tsv, monoploid_seqs) {

    def short_ids_head = short_ids_tsv.text.tokenize('\n')[0]

    if (short_ids_head == "IDs have acceptable length and character. No change required.") {
        return [ "${id}.mapped.monoploid.seqs.txt" ] + monoploid_seqs.text.tokenize('\n')
    }

    def orig_to_new_ids = [:]
    short_ids_tsv.text.eachLine { line ->
        def (original_id, renamed_id) = line.tokenize('\t')
        orig_to_new_ids[original_id] = renamed_id
    }

    def mapped_ids = []
    monoploid_seqs.text.eachLine { original_id ->
        if (!orig_to_new_ids[original_id]) {
            error "Faild to find $original_id in ${short_ids_tsv}" +
            "\nThe short_ids_tsv file is malformed!"
        }

        mapped_ids.add(orig_to_new_ids[original_id])
    }

    return [ "${id}.mapped.monoploid.seqs.txt" ] + mapped_ids
}


def restore_gff_ids(id, file_gff, short_ids_tsv) {
    def new_to_orig_ids = [:]

    short_ids_tsv.text.eachLine { line ->
        def columns = line.tokenize('\t')

        def orig_id = columns[0]
        def new_id  = columns[1]
        new_to_orig_ids[new_id] = orig_id
    }

    def restored_gff = []

    file_gff.text.eachLine { line ->
        if ( line.startsWith('##') ) {
            restored_gff.add(line)
            return
        }

        def columns = line.tokenize('\t')
        def new_id = columns[0]
        def orig_id = new_to_orig_ids[new_id]

        restored_gff.add ( ( [ orig_id ] + columns[1..-1] ).join('\t') )
    }

    [ "${id}.restored.ids.gff3", restored_gff.join('\n') ]

}


def shorten_ids(meta, fasta) {

    def ids_descs = fasta.splitFasta( record: [ id: true, desc: true ] )
    def ids = ids_descs.collect { it.id }

    if ( ids_descs.every { it.id.size() <= 13 && it.id ==~ /\w+/ && it.desc == null } ) {
        return [ "${meta.id}.short.ids.tsv", 'IDs have acceptable length and character. No change required.' ]
    }

    def shortened_ids = ids.collect { id -> ( id.size() <= 13 && id ==~ /\w+/ ) ? [ id, id ] : [ id, "Ctg${id.md5()[0..<10]}" ] }

    if ( shortened_ids.size() != shortened_ids.collect { it[1] }.unique().size() ) {
        error "Failed to create a unique set of Fasta IDs for sample ${meta.id}. Manual shortening of the Fasta IDs is required."
    } // The probability of this condition is very low but it is here to cover the edge cases
    // It assumes that the input fasta is valid and all the sequences have unique IDs

    [ "${meta.id}.short.ids.tsv", shortened_ids.collect { line -> line.join('\t') }.join('\n') ]
}
