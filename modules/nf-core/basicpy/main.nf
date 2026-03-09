process BASICPY {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/labsyspharm/basicpy-docker-mcmicro:1.2.0-patch5"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*-dfp.ome.tif"), path("*-ffp.ome.tif"), emit: profiles
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.2.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    /opt/main.py -i $image -o . --output-flatfield $prefix --output-darkfield $prefix $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basicpy: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.2.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    touch ${prefix}-dfp.ome.tif
    touch ${prefix}-ffp.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basicpy: $VERSION
    END_VERSIONS
    """
}
