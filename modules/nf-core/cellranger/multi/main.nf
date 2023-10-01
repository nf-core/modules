process CELLRANGER_MULTI {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:7.1.0"

    input:
    val meta
    tuple val(meta_gex)        , path (gex_fastqs   , stageAs: "fastqs/gex/*")
    tuple val(meta_vdj)        , path (vdj_fastqs   , stageAs: "fastqs/vdj/*")
    tuple val(meta_ab)         , path (ab_fastqs    , stageAs: "fastqs/ab/*")
    tuple val(meta_beam)       , path (beam_fastqs  , stageAs: "fastqs/beam/*")
    tuple val(meta_cmo)        , path (cmo_fastqs   , stageAs: "fastqs/cmo/*")
    tuple val(meta_crispr)     , path (crispr_fastqs, stageAs: "fastqs/crispr/*")
    path gex_reference         , stageAs: "references/gex/*"
    path gex_frna_probeset     , stageAs: "references/gex/probeset/*"
    path gex_targetpanel       , stageAs: "references/gex/targetpanel/*"
    path vdj_reference         , stageAs: "references/vdj/*"
    path vdj_primer_index      , stageAs: "references/vdj/primers/*"
    path fb_reference          , stageAs: "references/fb/*"
    path beam_panel            , stageAs: "references/beam/panel/*"
    path cmo_reference         , stageAs: "references/cmo/*"
    path cmo_barcodes          , stageAs: "references/cmo/barcodes/*"
    path cmo_barcode_assignment, stageAs: "references/cmo/sample_barcode_assignment/*"
    path frna_sampleinfo       , stageAs: "references/frna/*"

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
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // if references + FASTQ are empty, then don't run corresponding analyses
    // get names of references, if they exist
    // empty reference channels stage as "references"
    // empty FASTQ channels stage as "fastqs"
    // empty files stage as the file name, we check against 'EMPTY'
    def gex_reference_name     = gex_reference.getName() != 'references'     ? gex_reference.getName()          : ''
    def gex_frna_probeset_name = gex_frna_probeset.getBaseName() != 'EMPTY'  ? gex_frna_probeset.getName()      : ''
    def gex_targetpanel_name   = gex_targetpanel.getBaseName() != 'EMPTY'    ? gex_targetpanel.getName()        : ''
    def fb_reference_name      = fb_reference.getBaseName() != 'EMPTY'       ? fb_reference.getName()           : ''
    def vdj_reference_name     = vdj_reference.getName() != 'references'     ? vdj_reference.getName()          : ''
    def cmo_reference_name     = cmo_reference.getName() != 'EMPTY'          ? cmo_reference.getName()          : ''
    def cmo_sample_assignment  = cmo_barcode_assignment.getName() != 'EMPTY' ? cmo_barcode_assignment.getName() : ''

    def include_gex  = gex_fastqs.first().getName() != 'fastqs' && gex_reference ? '[gene-expression]'     : ''
    def include_vdj  = vdj_fastqs.first().getName() != 'fastqs' && vdj_reference ? '[vdj]'                 : ''
    def include_beam = beam_fastqs.first().getName() != 'fastqs' && beam_panel   ? '[antigen-specificity]' : ''
    def include_cmo  = cmo_fastqs.first().getName() != 'fastqs' && cmo_barcodes  ? '[samples]'             : ''
    def include_fb   = fb_reference.first().getName() != 'references'            ? '[feature]'             : ''
    def include_frna = gex_frna_probeset_name && frna_sampleinfo                 ? '[samples]'             : ''

    def gex_reference_path = include_gex ? "reference,\$PWD/${gex_reference_name}" : ''
    def fb_reference_path  = include_fb  ? "reference,\$PWD/${fb_reference_name}"  : ''
    def vdj_reference_path = include_vdj ? "reference,\$PWD/${vdj_reference_name}" : ''

    // targeted GEX panel goes under GEX section, not its own
    def target_panel = gex_targetpanel_name != '' ? "target-panel,\$PWD/$gex_targetpanel_name" : ''

    // fixed RNA reference (not sample info!) also goes under GEX section
    def frna_probeset = include_frna && gex_frna_probeset_name != '' ? "probe-set,\$PWD/references/gex/probeset/$gex_frna_probeset_name" : ''

    // VDJ inner primer set
    def primer_index = vdj_primer_index.getBaseName() != 'EMPTY' ? "inner-enrichment-primers,\$PWD/references/primers/${vdj_primer_index.getName()}" : ''

    // pull CSV text from these reference panels
    // these references get appended directly to config file
    def beam_csv_text  = include_beam && beam_panel.size() > 0      ? beam_panel.text      : ''
    def cmo_csv_text   = include_cmo  && cmo_barcodes.size() > 0    ? cmo_barcodes         : ''
    def frna_csv_text  = include_frna && frna_sampleinfo.size() > 0 ? frna_sampleinfo.text : ''

    // the feature barcodes section get options for either CRISPR or antibody capture assays
    def fb_options     = meta_ab?.options ? meta_ab.options : (meta_crispr?.options ? meta_crispr.options : [] )

    // collect options for each section
    // these are pulled from the meta maps
    def gex_options_use    = include_gex && meta_gex?.options   ? 'true' : null
    def vdj_options_use    = include_vdj && meta_vdj?.options   ? 'true' : null
    def ab_options_use     = include_fb && meta_ab?.options     ? 'true' : null
    def beam_options_use   = include_beam && meta_beam?.options ? 'true' : null
    def cmo_options_use    = include_cmo && meta_cmo?.options   ? 'true' : null
    def crispr_options_use = include_fb && meta_crispr?.options ? 'true' : null
    def fb_options_use     = include_fb && fb_options?.options  ? 'true' : null

    def gex_options_filter_probes = gex_options_use && meta_gex.options.containsKey("filter-probes") ? "filter-probes,${meta_gex.options["filter-probes"]}" : ''
    def gex_options_r1_length     = gex_options_use && meta_gex.options.containsKey("r1-length")     ? "r1-length,${meta_gex.options["r1-length"]}"         : ''
    def gex_options_r2_length     = gex_options_use && meta_gex.options.containsKey("r2-length")     ? "r2-length,${meta_gex.options["r2-length"]}"         : ''
    def gex_options_chemistry     = gex_options_use && meta_gex.options.containsKey("chemistry")     ? "chemistry,${meta_gex.options["chemistry"]}"         : ''
    def gex_options_expect_cells  = gex_options_use && meta_gex.options.containsKey("expect-cells")  ? "expect-cells,${meta_gex.options["expect-cells"]}"   : ''
    def gex_options_force_cells   = gex_options_use && meta_gex.options.containsKey("force-cells")   ? "force-cells,${meta_gex.options["force-cells"]}"     : ''
    def gex_options_no_secondary  = gex_options_use && meta_gex.options.containsKey("no-secondary")  ? "no-secondary,${meta_gex.options["no-secondary"]}"   : ''
    def gex_options_no_bam        = gex_options_use && meta_gex.options.containsKey("no-bam")        ? "no-bam,${meta_gex.options["no-bam"]}"               : ''
    def gex_options_no_target_umi_filter = gex_options_use && meta_gex.options.containsKey("no-target-umi-filter") ? "no-target-umi-filter,${meta_gex.options["no-target-umi-filter"]}" : ''
    def gex_options_include_introns      = gex_options_use && meta_gex.options.containsKey("include-introns")      ? "include-introns,${meta_gex.options["include-introns"]}"           : ''
    def gex_options_check_library_compatibility = gex_options_use && meta_gex.options.containsKey("check-library-compatibility") ? "check-library-compatibility,${meta_gex.options["check-library-compatibility"]}" : ''

    def cmo_reference_path = cmo_options_use && cmo_reference_name    ? "cmo-set,\$PWD/${cmo_reference_name}"                      : ''
    def cmo_barcode_path   = cmo_options_use && cmo_sample_assignment ? "barcode-sample-assignment,\$PWD/${cmo_sample_assignment}" : ''
    def cmo_options_min_assignment_confidence = cmo_options_use && meta_cmo.options.containsKey("min-assignment-confidence") ? "min-assignment-confidence,${meta_cmo.options["min-assignment-confidence"]}" : ''

    def vdj_options_r1_length = vdj_options_use && meta_vdj.options.containsKey("r1-length") ? "r1-length,${meta_vdj.options["r1-length"]}" : ''
    def vdj_options_r2_length = vdj_options_use && meta_vdj.options.containsKey("r2-length") ? "r2-length,${meta_vdj.options["r2-length"]}" : ''

    def fb_options_r1_length = fb_options_use && meta_fb.options.containsKey("r1-length") ? "r1-length,${meta_fb.options["r1-length"]}" : ''
    def fb_options_r2_length = fb_options_use && meta_fb.options.containsKey("r2-length") ? "r2-length,${meta_fb.options["r2-length"]}" : ''

    // point config to FASTQs
    def fastq_gex      = include_gex                      ? "${meta_gex.id},\$PWD/fastqs/gex/,,Gene Expression"           : ''
    def fastq_vdj      = include_vdj                      ? "${meta_vdj.id},\$PWD/fastqs/vdj,,VDJ"                        : ''
    def fastq_antibody = include_fb && ab_options_use     ? "${meta_ab.id},\$PWD/fastqs/ab,,Antibody Capture"             : ''
    def fastq_beam     = include_beam                     ? "${meta_beam.id},\$PWD/fastqs/beam,,Antigen Capture"          : ''
    def fastq_crispr   = include_fb && crispr_options_use ? "${meta_crispr.id},\$PWD/fastqs/crispr,,CRISPR Guide Capture" : ''
    def fastq_cmo      = include_cmo                      ? "${meta_cmo.id},\$PWD/fastqs/cmo,,Multiplexing Capture"       : ''

    // name the config file
    def config = "cellranger_multi_config.csv"

    """
    cat <<-CONFIG > $config
        $include_gex
        $gex_reference_path
        $frna_probeset
        $gex_options_filter_probes
        $gex_options_r1_length
        $gex_options_r2_length
        $gex_options_chemistry
        $gex_options_expect_cells
        $gex_options_force_cells
        $gex_options_no_secondary
        $gex_options_no_bam
        $gex_options_check_library_compatibility
        $target_panel
        $gex_options_no_target_umi_filter
        $gex_options_include_introns
        $cmo_options_min_assignment_confidence
        $cmo_reference_path
        $cmo_barcode_path

        $include_fb
        $fb_reference_path
        $fb_options_r1_length
        $fb_options_r2_length

        $include_vdj
        $vdj_reference_path
        $primer_index
        $vdj_options_r1_length
        $vdj_options_r2_length

        [libraries]
        fastq_id,fastqs,lanes,feature_types
        $fastq_gex
        $fastq_vdj
        $fastq_antibody
        $fastq_beam
        $fastq_crispr
        $fastq_cmo
    CONFIG

    if [[ "$include_cmo" ]]; then echo "$include_cmo" >> $config; fi
    if [[ "$include_cmo" ]]; then cat $cmo_barcodes >> $config; fi
    if [[ "$include_beam" ]]; then echo "$include_beam" >> $config; fi
    if [[ "$include_beam" ]]; then cat "$beam_csv_text" >> $config; fi
    if [[ "$include_frna" ]]; then echo "$include_frna" >> $config; fi
    if [[ "$include_frna" ]]; then cat "$frna_csv_text" >> $config; fi

    grep -v -e '^[[:space:]]*\$' $config > tmp.txt && mv tmp.txt $config # remove blank lines from config, only for aesthetics

    cellranger \\
        multi \\
        --id='${prefix}' \\
        --csv=$config \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
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
