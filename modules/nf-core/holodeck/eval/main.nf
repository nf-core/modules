process HOLODECK_EVAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fec611479d4b632095f4173ab61e9363051bd5645386a7873ce378846461cae9/data'
        : 'community.wave.seqera.io/library/holodeck:0.3.0--f29a8a9e2f667299'}"

    input:
    tuple val(meta), path(bam), path(truth)

    output:
    tuple val(meta), path("*.eval.txt"), emit: eval
    tuple val("${task.process}"), val('holodeck'), eval("holodeck --version | sed 's/^holodeck //'"), emit: versions_holodeck, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def truth_arg = truth ? "--truth ${truth}" : ""
    """
    holodeck \\
        eval \\
        --mapped ${bam} \\
        ${truth_arg} \\
        --output ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eval.txt
    """
}
