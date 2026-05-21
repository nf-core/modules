process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta),  path(tab), path(tai), path(regions)

    output:
    tuple val(meta), path("*.{tbi,csi}"),         emit: index,     optional: true
    tuple val(meta), path("${prefix}.*gz"),       emit: extracted, optional: true
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'")   , topic: versions   , emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix          = task.ext.prefix ?: "${meta.id}"
    def tab_suffix  = tab.name.indexOf('.') >= 0 ? tab.name.substring(tab.name.indexOf('.')) : ''
    def regions_arg = regions ? "-R ${regions}" : ""
    def output_arg  = regions ? "| bgzip --threads ${task.cpus} > ${prefix}${tab_suffix}" : ""
    """
    tabix \\
        ${regions_arg} \\
        --threads $task.cpus \\
        $args \\
        $tab \\
        ${output_arg}

    """
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tab_suffix = tab.name.indexOf('.') >= 0 ? tab.name.substring(tab.name.indexOf('.')) : ''
    def ext = args.contains("-C ") || args.contains("--csi") ? "csi" : "tbi"
    def index     = regions ? "" : "touch ${tab}.${ext}"
    def extracted = regions ? "echo | gzip > ${prefix}${tab_suffix}" : ""
    """
    ${index}
    ${extracted}
    """
}
