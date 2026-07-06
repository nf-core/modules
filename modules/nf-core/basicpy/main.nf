process BASICPY {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/labsyspharm/basicpy-docker-mcmicro:1.2.0-patch5"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*-dfp.ome.tif"), path("*-ffp.ome.tif"), emit: profiles
    tuple val("${task.process}"), val('basicpy'), val("1.2.0"), emit: versions_basicpy, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    /opt/main.py -i $image -o . --output-flatfield $prefix --output-darkfield $prefix $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-dfp.ome.tif
    touch ${prefix}-ffp.ome.tif
    """
}
