process BASICPY {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/yfukai/basicpy-docker-mcmicro:0.2.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.tiff"), emit: fields
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    /opt/main.py -i $image -o . $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basicpy:: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }
    """
    touch ${prefix}.-dfp.tiff
    touch ${prefix}.-dfp.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basicpy:: $VERSION
    END_VERSIONS
    """
}
