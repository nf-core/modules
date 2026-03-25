process GEMMI_VALIDATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/993d88cc40cbb0e8a42a8d173f7efdae2feb0924077d41c974c11459d49b6e5b/data':
        'community.wave.seqera.io/library/python_pip_gemmi-program:6276064ce54c78fb' }"

    input:
    tuple val(meta), path(cif)

    output:
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('gemmi'), eval("gemmi --version"), topic: versions, emit: versions_gemmi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gemmi \\
        validate \\
            $args \\
            ${cif} > ${prefix}.log 2>&1
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch ${prefix}.log
    """
}
