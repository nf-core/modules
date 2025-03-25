process CELLRANGER_MULTI {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:8.0.0"

    input:
    val meta
    tuple val(meta_gex)        , path (gex_fastqs   , stageAs: "fastqs/gex/fastq_???/*")
    tuple val(meta_vdj)        , path (vdj_fastqs   , stageAs: "fastqs/vdj/fastq_???/*")
    tuple val(meta_ab)         , path (ab_fastqs    , stageAs: "fastqs/ab/fastq_???/*")
    tuple val(meta_beam)       , path (beam_fastqs  , stageAs: "fastqs/beam/fastq_???/*")
    tuple val(meta_cmo)        , path (cmo_fastqs   , stageAs: "fastqs/cmo/fastq_???/*")
    tuple val(meta_crispr)     , path (crispr_fastqs, stageAs: "fastqs/crispr/fastq_???/*")
    path gex_reference         , stageAs: "references/gex/*"
    path gex_frna_probeset     , stageAs: "references/gex/probeset/*"
    path gex_targetpanel       , stageAs: "references/gex/targetpanel/*"
    path vdj_reference         , stageAs: "references/vdj/*"
    path vdj_primer_index      , stageAs: "references/vdj/primers/*"
    path fb_reference          , stageAs: "references/fb/*"
    path beam_antigen_panel    , stageAs: "references/beam/panel/antigens/*"
    path beam_control_panel    , stageAs: "references/beam/panel/controls/*"
    path cmo_reference         , stageAs: "references/cmo/*"
    path cmo_barcodes          , stageAs: "references/cmo/barcodes/*"
    path cmo_barcode_assignment, stageAs: "references/cmo/sample_barcode_assignment/*"
    path frna_sampleinfo       , stageAs: "references/frna/*"
    val skip_renaming

    output:
    tuple val(meta), path("cellranger_multi_config.csv"), emit: config
    tuple val(meta), path("**/outs/**")                 , emit: outs
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MULTI module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // if references + FASTQ are empty, then don't run corresponding analyses
    // get names of references, if they exist
    // empty reference channels (all under references/) can stage as "[]" when skipped by the workflow
    // empty FASTQ channels stage as "fastqs"
    gex_reference_name      = gex_reference          ? gex_reference.getName()          : ''
    gex_frna_probeset_name  = gex_frna_probeset      ? gex_frna_probeset.getName()      : ''
    gex_targetpanel_name    = gex_targetpanel        ? gex_targetpanel.getName()        : ''
    fb_reference_name       = fb_reference           ? fb_reference.getName()           : ''
    vdj_reference_name      = vdj_reference          ? vdj_reference.getName()          : ''
    cmo_reference_name      = cmo_reference          ? cmo_reference.getName()          : ''
    cmo_sample_assignment   = cmo_barcode_assignment ? cmo_barcode_assignment.getName() : ''
    beam_antigen_panel_name = beam_antigen_panel     ? beam_antigen_panel.getName()     : ''

    include_gex  = gex_fastqs.first().getName() != 'fastqs' && gex_reference           ? '[gene-expression]'     : ''
    include_vdj  = vdj_fastqs.first().getName() != 'fastqs' && vdj_reference           ? '[vdj]'                 : ''
    include_beam = beam_fastqs.first().getName() != 'fastqs' && beam_control_panel     ? '[antigen-specificity]' : ''
    include_cmo  = cmo_fastqs.first().getName() != 'fastqs' && cmo_barcodes            ? '[samples]'             : ''
    include_fb = (ab_fastqs.first().getName() != 'fastqs' || crispr_fastqs.first().getName() != 'fastqs') && fb_reference ? '[feature]' : ''
    include_frna = gex_frna_probeset_name && frna_sampleinfo                           ? '[samples]'             : ''

    gex_reference_path = include_gex ? "reference,./${gex_reference_name}" : ''
    fb_reference_path  = include_fb  ? "reference,./${fb_reference_name}"  : ''
    vdj_reference_path = include_vdj ? "reference,./${vdj_reference_name}" : ''

    // targeted GEX panel goes under GEX section, not its own
    target_panel = gex_targetpanel_name != '' ? "target-panel,./$gex_targetpanel_name" : ''

    // fixed RNA reference (not sample info!) also goes under GEX section
    frna_probeset = include_frna && gex_frna_probeset_name != '' ? "probe-set,./$gex_frna_probeset_name" : ''

    // VDJ inner primer set
    primer_index = vdj_primer_index ? "inner-enrichment-primers,./references/primers/${vdj_primer_index.getName()}" : ''

    // BEAM antigen list, remember that this is a Feature Barcode file
    beam_antigen_csv = include_beam && beam_antigen_panel_name != '' ? "reference,./$beam_antigen_panel_name" : ''

    // pull CSV text from these reference panels
    // these references get appended directly to config file
    beam_csv_text  = include_beam && beam_control_panel.size() > 0 ? beam_control_panel : ''
    cmo_csv_text   = include_cmo  && cmo_barcodes.size() > 0       ? cmo_barcodes       : ''
    frna_csv_text  = include_frna && frna_sampleinfo.size() > 0    ? frna_sampleinfo    : ''

    // the feature barcodes section get options for either CRISPR or antibody capture assays
    fb_options     = meta_ab?.options ? meta_ab.options : (meta_crispr?.options ? meta_crispr.options : [])

    // collect options for each section
    // these are pulled from the meta maps
    gex_options_use    = include_gex && meta_gex?.options   ? 'true' : null
    vdj_options_use    = include_vdj && meta_vdj?.options   ? 'true' : null
    ab_options_use     = include_fb && meta_ab?.options     ? 'true' : null
    beam_options_use   = include_beam && meta_beam?.options ? 'true' : null
    cmo_options_use    = include_cmo && meta_cmo?.options   ? 'true' : null
    crispr_options_use = include_fb && meta_crispr?.options ? 'true' : null
    fb_options_use     = include_fb && fb_options?.options  ? 'true' : null

    gex_options_filter_probes = gex_options_use && meta_gex.options.containsKey("filter-probes") ? "filter-probes,${meta_gex.options["filter-probes"]}" : ''
    gex_options_r1_length     = gex_options_use && meta_gex.options.containsKey("r1-length")     ? "r1-length,${meta_gex.options["r1-length"]}"         : ''
    gex_options_r2_length     = gex_options_use && meta_gex.options.containsKey("r2-length")     ? "r2-length,${meta_gex.options["r2-length"]}"         : ''
    gex_options_chemistry     = gex_options_use && meta_gex.options.containsKey("chemistry")     ? "chemistry,${meta_gex.options["chemistry"]}"         : ''
    gex_options_expect_cells  = gex_options_use && meta_gex.options.containsKey("expect-cells")  ? "expect-cells,${meta_gex.options["expect-cells"]}"   : ''
    gex_options_force_cells   = gex_options_use && meta_gex.options.containsKey("force-cells")   ? "force-cells,${meta_gex.options["force-cells"]}"     : ''
    gex_options_no_secondary  = gex_options_use && meta_gex.options.containsKey("no-secondary")  ? "no-secondary,${meta_gex.options["no-secondary"]}"   : ''
    gex_options_no_bam        = gex_options_use && meta_gex.options.containsKey("create-bam")    ? "create-bam,${meta_gex.options["create-bam"]}"           : ''
    gex_options_no_target_umi_filter = gex_options_use && meta_gex.options.containsKey("no-target-umi-filter") ? "no-target-umi-filter,${meta_gex.options["no-target-umi-filter"]}" : ''
    gex_options_include_introns      = gex_options_use && meta_gex.options.containsKey("include-introns")      ? "include-introns,${meta_gex.options["include-introns"]}"           : ''
    gex_options_check_library_compatibility = gex_options_use && meta_gex.options.containsKey("check-library-compatibility") ? "check-library-compatibility,${meta_gex.options["check-library-compatibility"]}" : ''

    cmo_reference_path = cmo_options_use && cmo_reference_name    ? "cmo-set,./${cmo_reference_name}"                      : ''
    cmo_barcode_path   = cmo_options_use && cmo_sample_assignment ? "barcode-sample-assignment,./${cmo_sample_assignment}" : ''
    cmo_options_min_assignment_confidence = cmo_options_use && meta_cmo.options.containsKey("min-assignment-confidence") ? "min-assignment-confidence,${meta_cmo.options["min-assignment-confidence"]}" : ''

    vdj_options_r1_length = vdj_options_use && meta_vdj.options.containsKey("r1-length") ? "r1-length,${meta_vdj.options["r1-length"]}" : ''
    vdj_options_r2_length = vdj_options_use && meta_vdj.options.containsKey("r2-length") ? "r2-length,${meta_vdj.options["r2-length"]}" : ''

    fb_options_r1_length = fb_options_use && meta_fb.options.containsKey("r1-length") ? "r1-length,${meta_fb.options["r1-length"]}" : ''
    fb_options_r2_length = fb_options_use && meta_fb.options.containsKey("r2-length") ? "r2-length,${meta_fb.options["r2-length"]}" : ''

    // point config to FASTQs
    // After renaming it gets in 'fastq_all' folder
    fastq_gex      = include_gex                      ? "${meta_gex.id},./fastq_all/gex,,Gene Expression"            : ''
    fastq_vdj      = include_vdj                      ? "${meta_vdj.id},./fastq_all/vdj,,VDJ"                        : ''
    fastq_antibody = include_fb && ab_options_use     ? "${meta_ab.id},./fastq_all/ab,,Antibody Capture"             : ''
    fastq_beam     = include_beam                     ? "${meta_beam.id},./fastq_all/beam,,Antigen Capture"         : ''
    fastq_crispr   = include_fb && crispr_options_use ? "${meta_crispr.id},./fastq_all/crispr,,CRISPR Guide Capture" : ''
    fastq_cmo      = include_cmo                      ? "${meta_cmo.id},./fastq_all/cmo,,Multiplexing Capture"       : ''

    // name the config file
    config = "cellranger_multi_config.csv"
    template "cellranger_multi.py"

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    touch ${prefix}/outs/fake_file.txt
    echo -n "" >> ${prefix}/outs/fake_file.txt
    touch cellranger_multi_config.csv
    echo -n "" >> cellranger_multi_config.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
