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
    tuple val("${task.process}"), val('simpleaf'),   eval("simpleaf --version | sed 's/simpleaf //'"),                          topic: versions, emit: versions_simpleaf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def (forward, reverse) = reads.collate(2).transpose()
    def reads1 = forward.join(',')
    def reads2 = reverse.join(',')

    def index_arg          = index          ? "--index ${index}"                     : ''
    def probe_set_arg      = probe_set      ? "--probe-set ${probe_set}"             : ''
    def sample_bc_list_arg = sample_bc_list ? "--sample-bc-list ${sample_bc_list}"   : ''
    def cell_bc_list_arg   = cell_bc_list   ? "--cell-bc-list ${cell_bc_list}"       : ''
    def t2g_map_arg        = t2g_map        ? "--t2g-map ${t2g_map}"                 : ''

    meta = meta2 + meta3 + meta

    """
    export ALEVIN_FRY_HOME=.
    simpleaf set-paths

    simpleaf multiplex-quant \\
        --chemistry ${chemistry} \\
        --reads1 ${reads1} \\
        --reads2 ${reads2} \\
        ${index_arg} \\
        ${probe_set_arg} \\
        ${sample_bc_list_arg} \\
        ${cell_bc_list_arg} \\
        ${t2g_map_arg} \\
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
