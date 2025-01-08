process BAYSOR_SEGFREE {
    tag '$meta.id'
    label 'process_high_memory'
    container "docker.io/segonzal/baysor:0.7.1"

    input:
    tuple val(meta), path(transcripts_csv)
    path(config_toml)

    output:
    tuple val(meta), path("ncvs.loom"), emit: ncvs
    path("ncvs_segfree_log.log"),       emit: ncvs_log

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BAYSOR_SEGFREE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.7.1"

    """
    baysor segfree ${transcripts_csv} -c ${config_toml} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BAYSOR_SEGFREE module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.7.1"

    """
    touch ncvs.loom
    touch ncvs_segfree_log.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor: $VERSION
    END_VERSIONS
    """
}
