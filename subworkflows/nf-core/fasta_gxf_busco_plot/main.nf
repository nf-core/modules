include { BUSCO_BUSCO as BUSCO_ASSEMBLY         } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as PLOT_ASSEMBLY   } from '../../../modules/nf-core/busco/generateplot/main'
include { GFFREAD as EXTRACT_PROTEINS           } from '../../../modules/nf-core/gffread/main'
include { BUSCO_BUSCO as BUSCO_ANNOTATION       } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as PLOT_ANNOTATION } from '../../../modules/nf-core/busco/generateplot/main'

workflow FASTA_GXF_BUSCO_PLOT {

    take:
    ch_fasta                                    // channel: [ val(meta), fasta ]
    ch_gxf                                      // channel: [ val(meta2), gxf ]; gxf ~ gff | gff3 | gtf
                                                //
                                                // meta and meta2 should have same id

    val_mode                                    // val(mode); BUSCO mode to apply to ch_fasta
                                                // - genome, for genome assemblies (DNA)
                                                // - transcriptome, for transcriptome assemblies (DNA)
                                                // - proteins, for annotated gene sets (protein)
                                                //
                                                // If mode is genome, annotations from ch_gxf are evaluated with
                                                // mode proteins, otherwise, evaluation of the annotations is skipped
                                                //
    val_lineages                                // [ val(lineage) ]
    val_busco_lineages_path                     // val(path); Optional; Set to [] if not needed
    val_busco_config                            // val(path); Optional; Set to [] if not needed

    main:
    ch_versions                                 = Channel.empty()
    ch_db_path                                  = val_busco_lineages_path
                                                ? Channel.of(file(val_busco_lineages_path, checkIfExists: true))
                                                : Channel.of( [ [] ] )
    ch_config_path                              = val_busco_config
                                                ? Channel.of(file(val_busco_config, checkIfExists: true))
                                                : Channel.of( [ [] ] )

    // MODULE: BUSCO_BUSCO as BUSCO_ASSEMBLY
    ch_busco_assembly_inputs                    = ch_fasta
                                                | combine(
                                                    Channel.of(val_mode)
                                                )
                                                | combine(
                                                    Channel.fromList(val_lineages)
                                                )
                                                | map { meta, fasta, mode, lineage ->
                                                    [
                                                        meta + [ mode:mode, lineage:lineage ],
                                                        fasta, mode, lineage
                                                    ]
                                                }
                                                | combine(
                                                    ch_db_path
                                                )
                                                | combine(
                                                    ch_config_path
                                                )

    BUSCO_ASSEMBLY(
        ch_busco_assembly_inputs.map { meta,  fasta,  _mode, _lineage, _db, _config -> [ meta, fasta ] },
        ch_busco_assembly_inputs.map { _meta, _fasta, mode,  _lineage, _db, _config -> mode },
        ch_busco_assembly_inputs.map { _meta, _fasta, _mode, lineage,  _db, _config -> lineage },
        ch_busco_assembly_inputs.map { _meta, _fasta, _mode, _lineage, db,  _config -> db },
        ch_busco_assembly_inputs.map { _meta, _fasta, _mode, _lineage, _db, config  -> config }
    )

    ch_assembly_batch_summary                   = BUSCO_ASSEMBLY.out.batch_summary
    ch_assembly_short_summaries_txt             = BUSCO_ASSEMBLY.out.short_summaries_txt
    ch_assembly_short_summaries_json            = BUSCO_ASSEMBLY.out.short_summaries_json
    ch_assembly_full_table                      = BUSCO_ASSEMBLY.out.full_table
    ch_versions                                 = ch_versions.mix(BUSCO_ASSEMBLY.out.versions.first())

    // MODULE: BUSCO_GENERATEPLOT as PLOT_ASSEMBLY
    ch_assembly_plot_summary                    = ch_assembly_short_summaries_txt
                                                | map { meta, txt ->
                                                    def lineage_name = meta.lineage - '_odb'
                                                    [
                                                        "short_summary.specific.${meta.lineage}.${meta.id}_${lineage_name}.txt",
                                                        txt.text
                                                    ]
                                                }
                                                | collectFile

    PLOT_ASSEMBLY( ch_assembly_plot_summary.collect() )

    ch_assembly_png                             = PLOT_ASSEMBLY.out.png
    ch_versions                                 = ch_versions.mix(PLOT_ASSEMBLY.out.versions)

    // MODULE: GFFREAD as EXTRACT_PROTEINS
    ch_gffread_inputs                           = val_mode !in [ 'geno', 'genome' ]
                                                ? Channel.empty()
                                                : ch_fasta
                                                | map { meta, fasta -> [ meta.id, meta, fasta ] }
                                                | join(
                                                    ch_gxf.map { meta2, gxf -> [ meta2.id, gxf ] }
                                                    // Join with matching annotation
                                                    // to allow one annotations per fasta
                                                )
                                                | map { _id, meta, fasta, gxf -> [ meta, gxf, fasta ] }
    EXTRACT_PROTEINS(
        ch_gffread_inputs.map { meta,  gxf,  _fasta -> [ meta, gxf ] },
        ch_gffread_inputs.map { _meta, _gxf, fasta  -> fasta }
    )

    ch_proteins                                 = EXTRACT_PROTEINS.out.gffread_fasta
    ch_versions                                 = ch_versions.mix(EXTRACT_PROTEINS.out.versions.first())

    // MODULE: BUSCO_BUSCO as BUSCO_ANNOTATION
    ch_busco_annotation_inputs                  = ch_proteins
                                                | combine(
                                                    Channel.of('proteins')
                                                )
                                                | combine(
                                                    Channel.fromList(val_lineages)
                                                )
                                                | map { meta, fasta, mode, lineage ->
                                                    [
                                                        meta + [ mode:mode, lineage:lineage ],
                                                        fasta, mode, lineage
                                                    ]
                                                }
                                                | combine(
                                                    ch_db_path
                                                )
                                                | combine(
                                                    ch_config_path
                                                )

    BUSCO_ANNOTATION(
        ch_busco_annotation_inputs.map { meta,  fasta,  _mode, _lineage, _db, _config -> [ meta, fasta ] },
        ch_busco_annotation_inputs.map { _meta, _fasta, mode,  _lineage, _db, _config -> mode },
        ch_busco_annotation_inputs.map { _meta, _fasta, _mode, lineage,  _db, _config -> lineage },
        ch_busco_annotation_inputs.map { _meta, _fasta, _mode, _lineage, db,  _config -> db },
        ch_busco_annotation_inputs.map { _meta, _fasta, _mode, _lineage, _db, config  -> config }
    )

    ch_annotation_batch_summary                 = BUSCO_ANNOTATION.out.batch_summary
    ch_annotation_short_summaries_txt           = BUSCO_ANNOTATION.out.short_summaries_txt
    ch_annotation_short_summaries_json          = BUSCO_ANNOTATION.out.short_summaries_json
    ch_annotation_full_table                    = BUSCO_ANNOTATION.out.full_table
    ch_versions                                 = ch_versions.mix(BUSCO_ANNOTATION.out.versions.first())

    // MODULE: BUSCO_GENERATEPLOT as PLOT_ANNOTATION
    ch_annotation_plot_summary                  = ch_annotation_short_summaries_txt
                                                | map { meta, txt ->
                                                    def lineage_name = meta.lineage - '_odb'
                                                    [
                                                        "short_summary.specific.${meta.lineage}.${meta.id}_${lineage_name}.proteins.txt",
                                                        txt.text
                                                    ]
                                                }
                                                | collectFile

    PLOT_ANNOTATION( ch_annotation_plot_summary.collect() )

    ch_annotation_png                           = PLOT_ANNOTATION.out.png
    ch_versions                                 = ch_versions.mix(PLOT_ANNOTATION.out.versions)


    emit:
    assembly_batch_summary                      = ch_assembly_batch_summary             // channel: [ meta3, txt ]; meta3 ~ meta + [ val(mode), val(lineage) ]
    assembly_short_summaries_txt                = ch_assembly_short_summaries_txt       // channel: [ meta3, txt ]
    assembly_short_summaries_json               = ch_assembly_short_summaries_json      // channel: [ meta3, json ]
    assembly_full_table                         = ch_assembly_full_table                // channel: [ meta3, tsv ]
    assembly_plot_summary_txt                   = ch_assembly_plot_summary              // channel: [ text ]
    assembly_png                                = ch_assembly_png                       // channel: [ png ]
    annotation_batch_summary                    = ch_annotation_batch_summary           // channel: [ meta3, txt ]
    annotation_short_summaries_txt              = ch_annotation_short_summaries_txt     // channel: [ meta3, txt ]
    annotation_short_summaries_json             = ch_annotation_short_summaries_json    // channel: [ meta3, json ]
    annotation_full_table                       = ch_annotation_full_table              // channel: [ meta3, tsv ]
    annotation_plot_summary_txt                 = ch_annotation_plot_summary            // channel: [ txt ]
    annotation_png                              = ch_annotation_png                     // channel: [ png ]
    versions                                    = ch_versions                           // channel: [ versions.yml ]
}
