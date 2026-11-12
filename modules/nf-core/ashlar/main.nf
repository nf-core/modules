// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ASHLAR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5f9a902bb541dc2abba28e73c4bac2ee14d28373abe2e14b444a5096f4919cb7/data'
        : 'community.wave.seqera.io/library/ashlar_tifffile:6466db6f351397ba'}"

    input:
    tuple val(meta), path(images, stageAs: 'image*/*')
    path(opt_dfp, stageAs: 'dfp*/*')
    path(opt_ffp, stageAs: 'ffp*/*')

    output:
    tuple val(meta), path("*.ome.tif"), emit: tif
    tuple val("${task.process}"), val('ashlar'), eval("ashlar --version | sed 's/^.*ashlar //'"), emit: versions_ashlar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args           ?: ''
    def prefix        = task.ext.prefix         ?: "${meta.id}"
    def dfp           = opt_dfp                 ? "--dfp ${opt_dfp}" : ""
    def ffp           = opt_ffp                 ? "--ffp ${opt_ffp}" : ""
    def num_files     = images instanceof List  ? images.size()      : 1
    def opt_dfp_size  = opt_dfp instanceof List ? opt_dfp.size()     : 1
    def opt_ffp_size  = opt_ffp instanceof List ? opt_ffp.size()     : 1
    def dfp_validated = (opt_dfp_size == 0 || opt_dfp_size == 1 || opt_dfp_size == num_files) ? true : false
    def ffp_validated = (opt_ffp_size == 0 || opt_ffp_size == 1 || opt_ffp_size == num_files) ? true : false

    if ( !dfp_validated ) { error "Please input only zero, one, or N dfp files, where N is the number of input images" }
    if ( !ffp_validated ) { error "Please input only zero, one, or N ffp files, where N is the number of input images" }

    """
    export JAVA_TOOL_OPTIONS='-XX:+PerfDisableSharedMem'
    export XDG_CACHE_HOME=.cache
    export MPLCONFIGDIR=.matplotlib

    ashlar \\
        -o ${prefix}.ome.tif \\
        ${images} \\
        ${args} \\
        ${dfp} \\
        ${ffp}

    sed -i -E 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' ${prefix}.ome.tif
    """

    stub:
    def prefix        = task.ext.prefix         ?: "${meta.id}"
    def num_files     = images instanceof List  ? images.size()  : 1
    def opt_dfp_size  = opt_dfp instanceof List ? opt_dfp.size() : 1
    def opt_ffp_size  = opt_ffp instanceof List ? opt_ffp.size() : 1
    def dfp_validated = (opt_dfp_size == 0 || opt_dfp_size == 1 || opt_dfp_size == num_files) ? true : false
    def ffp_validated = (opt_ffp_size == 0 || opt_ffp_size == 1 || opt_ffp_size == num_files) ? true : false

    if ( !dfp_validated ) { error "Please input only zero, one, or N dfp files, where N is the number of input images" }
    if ( !ffp_validated ) { error "Please input only zero, one, or N ffp files, where N is the number of input images" }

    """
    export JAVA_TOOL_OPTIONS='-XX:+PerfDisableSharedMem'

    touch ${prefix}.ome.tif
    """
}
