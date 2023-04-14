process BASICPY {
    tag "$meta.id"
    label 'process_single'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Basicpy module does not support Conda. Please use Docker / Singularity instead."
    }

    container "yfukai/basicpy-docker-mcmicro:0.1.2"

    input:
    tuple val(meta), path(image)
    val(cpu_gpu)

    output:
    tuple val(meta), path("*.tiff"), emit: fields
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.1.2" // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    /opt/main.py $cpu_gpu $image . $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basicpy:: $VERSION
    END_VERSIONS
    """
}
