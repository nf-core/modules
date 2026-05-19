process SIMPLEAF_MULTIPLEXQUANT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/aa/aaba033a0179fd6ccc20c677f9df1fac5d8eac2dbd1bed73c4fa9f7adb65d963/data':
        'community.wave.seqera.io/library/simpleaf:0.25.0--b9f96d8b71a01864' }"

    input:
    //
    // Input reads are expected as: [ meta, chemistry_preset, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Reads are split into R1/R2 pairs and joined with commas before being passed to simpleaf.
    //
    tuple val(meta),  val(chemistry), path(reads)                                                           // chemistry preset and reads
    tuple val(meta2), path(index, stageAs: 'index/*'), path(t2g_map)                                        // optional pre-built piscem probe index and t2g map
    tuple val(meta3), path(probe_set), path(sample_bc_list), path(cell_bc_list)                             // optional probe set / sample-BC TSV / cell-BC whitelist overrides
    val resolution                                                                                          // UMI resolution (cr-like, cr-like-em, parsimony, ...)

    output:
    tuple val(meta), path("${prefix}/af_map")                       , emit: map
    tuple val(meta), path("${prefix}/af_quant")                     , emit: quant
    tuple val(meta), path("${prefix}/af_quant/alevin/quants.h5ad")  , emit: h5ad,        optional: true
    tuple val(meta), path("${prefix}/probe_t2g.tsv")                , emit: t2g,         optional: true
    tuple val(meta), path("${prefix}/probe_index/index")            , emit: probe_index, optional: true
    tuple val("${task.process}"), val('alevin-fry'), eval("alevin-fry --version | sed 's/alevin-fry //'"),                    topic: versions, emit: versions_alevin_fry
    tuple val("${task.process}"), val('piscem'),     eval("piscem --version | sed 's/piscem //'"),                            topic: versions, emit: versions_piscem
    tuple val("${task.process}"), val('simpleaf'),   eval("ALEVIN_FRY_HOME=. simpleaf --version | sed 's/simpleaf //'"),      topic: versions, emit: versions_simpleaf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def mapping_args   = mappingArgs(chemistry, reads)
    def reference_args = referenceArgs(index, probe_set, sample_bc_list, cell_bc_list, t2g_map)

    meta = meta2 + meta3 + meta

    """
    export ALEVIN_FRY_HOME=.
    simpleaf set-paths

    # run simpleaf multiplex-quant
    simpleaf multiplex-quant \\
        ${mapping_args} \\
        ${reference_args} \\
        --resolution ${resolution} \\
        --output ${prefix} \\
        --threads ${task.cpus} \\
        --anndata-out \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export ALEVIN_FRY_HOME=.

    mkdir -p ${prefix}/af_map
    mkdir -p ${prefix}/af_quant/alevin

    touch ${prefix}/af_map/map.rad
    touch ${prefix}/af_map/map_info.json
    touch ${prefix}/af_quant/quant.json
    touch ${prefix}/af_quant/generate_permit_list.json
    touch ${prefix}/af_quant/alevin/quants_mat.mtx
    touch ${prefix}/af_quant/alevin/quants_mat_rows.txt
    touch ${prefix}/af_quant/alevin/quants_mat_cols.txt
    touch ${prefix}/af_quant/alevin/quants.h5ad
    touch ${prefix}/probe_t2g.tsv
    """
}

// `simpleaf multiplex-quant` requires both reads and a chemistry preset (or, with extra
// ext.args, a --geometry override + --cell-bc-list). Only the mainstream case is enforced
// here; non-default geometries can still be set via ext.args.
def mappingArgs(chemistry, reads) {
    if (!reads)     error "Missing read files; could not proceed."
    if (!chemistry) error "Missing chemistry; could not proceed."

    def (forward, reverse) = reads.collate(2).transpose()
    return """--chemistry ${chemistry} \\
        --reads1 ${forward.join(',')} \\
        --reads2 ${reverse.join(',')}"""
}

// Build optional reference-override flags. With none of these set, simpleaf auto-downloads
// a probe set + sample BC TSV based on the chemistry preset and (if also provided in ext.args)
// `--organism`. Any combination of overrides is allowed.
def referenceArgs(index, probe_set, sample_bc_list, cell_bc_list, t2g_map) {
    def parts = []
    if (index)          parts << "--index ${index}"
    if (probe_set)      parts << "--probe-set ${probe_set}"
    if (sample_bc_list) parts << "--sample-bc-list ${sample_bc_list}"
    if (cell_bc_list)   parts << "--cell-bc-list ${cell_bc_list}"
    if (t2g_map)        parts << "--t2g-map ${t2g_map}"
    return parts.join(' \\\n        ')
}
