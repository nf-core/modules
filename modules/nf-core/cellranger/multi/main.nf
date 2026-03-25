process CELLRANGER_MULTI {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/cellranger:10.0.0"

    input:
    val meta
    tuple val(meta2)           , path (gex_fastqs   , stageAs: "fastqs/gex/fastq_???/*")
    tuple val(meta3)           , path (vdj_fastqs   , stageAs: "fastqs/vdj/fastq_???/*")
    tuple val(meta4)           , path (ab_fastqs    , stageAs: "fastqs/ab/fastq_???/*")
    tuple val(meta5)           , path (beam_fastqs  , stageAs: "fastqs/beam/fastq_???/*")
    tuple val(meta6)           , path (cmo_fastqs   , stageAs: "fastqs/cmo/fastq_???/*")
    tuple val(meta7)           , path (crispr_fastqs, stageAs: "fastqs/crispr/fastq_???/*")
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
    path ocm_barcodes          , stageAs: "references/ocm/barcodes/*"
    val skip_renaming

    output:
    tuple val(meta), path("cellranger_multi_config.csv"), emit: config
    tuple val(meta), path("**/outs/**")                 , emit: outs
    tuple val("${task.process}"), val('cellranger'), eval('cellranger --version | sed "s/.* //"'), emit: versions_cellranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MULTI module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    
    // Validate mutually exclusive barcode types
    if ([ocm_barcodes, cmo_barcodes, frna_sampleinfo].findAll().size() >= 2) {
        error "The ocm barcodes, cmo barcodes, and frna probes are mutually exclusive features. Please use only one per sample."
    }
    
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: meta.id

    // Determine which library types are present based on FASTQs and references
    def has_gex    = meta2 && gex_reference
    def has_vdj    = meta3 && vdj_reference
    def has_ab     = meta4 && fb_reference
    def has_beam   = meta5 && beam_control_panel
    def has_cmo    = meta6 && cmo_barcodes
    def has_crispr = meta7 && fb_reference
    def has_frna   = gex_frna_probeset && frna_sampleinfo
    def has_ocm    = ocm_barcodes

    // Build [gene-expression] section
    def gex_section = []
    if (has_gex) {
        gex_section << '[gene-expression]'
        gex_section << "reference,\$PWD/${gex_reference.name}"
        if (gex_frna_probeset) gex_section << "probe-set,\$PWD/${gex_frna_probeset.name}"
        
        // GEX options from meta
        def gex_has_opts = meta2?.options
        def chemistry = gex_has_opts?.containsKey("chemistry") ? meta2.options["chemistry"] : "auto"
        gex_section << "chemistry,${chemistry}"
        
        if (gex_has_opts?.containsKey("expect-cells")) {
            gex_section << "expect-cells,${meta2.options["expect-cells"]}"
        }
        
        def create_bam = gex_has_opts?.containsKey("create-bam") ? meta2.options["create-bam"] : "true"
        gex_section << "create-bam,${create_bam}"
        
        if (gex_targetpanel) {
            gex_section << "target-panel,\$PWD/${gex_targetpanel.name}"
        }
    }

    // Build [feature] section
    def fb_section = []
    if (has_ab || has_crispr) {
        fb_section << '[feature]'
        fb_section << 'reference,\$PWD/fb_reference_copy.csv'
    }

    // Build [vdj] section
    def vdj_section = []
    if (has_vdj) {
        vdj_section << '[vdj]'
        vdj_section << "reference,\$PWD/${vdj_reference.name}"
        if (vdj_primer_index) {
            vdj_section << "inner-enrichment-primers,\$PWD/${vdj_primer_index.name}"
        }
    }

    // Build [libraries] section
    def lib_section = ['[libraries]', 'fastq_id,fastqs,lanes,feature_types']
    if (has_gex) lib_section << "${meta2.id},\$PWD/fastq_all/gex,,Gene Expression"
    if (has_vdj) lib_section << "${meta3.id},\$PWD/fastq_all/vdj,,VDJ"
    if (has_ab) lib_section << "${meta4.id},\$PWD/fastq_all/ab,,Antibody Capture"
    if (has_beam) lib_section << "${meta5.id},\$PWD/fastq_all/beam,,Antigen Capture"
    if (has_crispr) lib_section << "${meta7.id},\$PWD/fastq_all/crispr,,CRISPR Guide Capture"
    if (has_cmo) lib_section << "${meta6.id},\$PWD/fastq_all/cmo,,Multiplexing Capture"

    // Build config content by combining all sections
    def config_lines = []
    config_lines.addAll(gex_section)
    config_lines.addAll(fb_section)
    config_lines.addAll(vdj_section)
    config_lines.addAll(lib_section)
    
    // Append sample sections if present
    if (has_cmo) {
        config_lines << '[samples]'
        config_lines << cmo_barcodes.text.trim()
    }
    if (has_frna) {
        config_lines << '[samples]'
        config_lines << frna_sampleinfo.text.trim()
    }
    if (has_ocm) {
        config_lines << '[samples]'
        config_lines << ocm_barcodes.text.trim()
    }
    def config_content = config_lines.findAll { line -> line }.join('\n    ')
    """
    #
    # Rename FASTQs to Cell Ranger naming convention
    # Maintains R1/R2 pairing order by processing directories sequentially
    # Output pattern: \${prefix}_S1_L001_R1_001.fastq.gz, L002_R1_001.fastq.gz, etc.
    #
    mkdir -p fastq_all/{gex,vdj,ab,beam,cmo,crispr}

    for modality in gex vdj ab beam cmo crispr; do
        lane=1
        while IFS= read -r -d '' r1_dir && IFS= read -r -d '' r2_dir; do
            r1=\$(find "\${r1_dir}" -maxdepth 1 -name "*_R1_*.fastq.gz" | head -1)
            r2=\$(find "\${r2_dir}" -maxdepth 1 -name "*_R2_*.fastq.gz" | head -1)
            
            [ -z "\${r1}" ] || [ -z "\${r2}" ] && continue
            
            ln -sf "\$(readlink -f "\${r1}")" "fastq_all/\${modality}/${prefix}_S1_L\$(printf %03d \${lane})_R1_001.fastq.gz"
            ln -sf "\$(readlink -f "\${r2}")" "fastq_all/\${modality}/${prefix}_S1_L\$(printf %03d \${lane})_R2_001.fastq.gz"
            lane=\$((lane + 1))
        done < <(find fastqs/\${modality} -maxdepth 1 -type d -name "fastq_*" | sort | xargs -n1 printf '%s\\0')
    done

    #
    # Copy fb_reference to avoid symlink corruption
    # Cell Ranger writes to this file during validation, which corrupts the symlinked original
    #
    if [ -n "${fb_reference}" ] && [ -f "${fb_reference}" ]; then
        cp "${fb_reference}" "fb_reference_copy.csv"
    fi

    cat > cellranger_multi_config.csv <<-CONFIG_EOF
    ${config_content}
    CONFIG_EOF

    cellranger multi \\
        --id=${prefix} \\
        --csv=cellranger_multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}/outs/"
    touch ${prefix}/outs/fake_file.txt
    echo -n "" >> ${prefix}/outs/fake_file.txt
    touch cellranger_multi_config.csv
    echo -n "" >> cellranger_multi_config.csv
    """
}
