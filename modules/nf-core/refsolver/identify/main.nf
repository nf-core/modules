process REFSOLVER_IDENTIFY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d900b9428fe55a93e8beeb3d01e594d683dc8526e26fc1105f33738ac6ad698/data'
        : 'community.wave.seqera.io/library/ref-solver:0.3.0--f7277dbf424ce3e5'}"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: report
    tuple val("${task.process}"), val('refsolver'), eval("ref-solver --version | sed 's/^ref-solver //'"), emit: versions_refsolver, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = args.contains("json") ? "json" : args.contains("tsv") ? "tsv" : "txt"
    """
    ref-solver \\
        identify \\
        ${args} \\
        ${alignment} \\
        > ${prefix}.${suffix}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = args.contains("json") ? "json" : args.contains("tsv") ? "tsv" : "txt"
    """
    touch ${prefix}.${suffix}
    """
}
