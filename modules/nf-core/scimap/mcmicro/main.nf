process SCIMAP_MCMICRO {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/labsyspharm/scimap:0.22.0"

    input:
    tuple val(meta), path(cellbyfeature)

    output:
    tuple val(meta), path("*.csv")      , emit: annotedDataCsv, optional:true
    tuple val(meta), path("*.h5ad")     , emit: annotedDataH5ad, optional:true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Scimap module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION='0.22.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping
    """
    scimap-mcmicro $cellbyfeature -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scimap: $VERSION
    END_VERSIONS
    """

}
