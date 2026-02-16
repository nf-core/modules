process WINDOWMASKER_MKCOUNTS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*.txt")  , emit: counts
    tuple val("${task.process}"), val('windowmasker'), eval("windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker. //; s/ .*\$//'"), topic: versions, emit: versions_windowmasker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    def memory  = 3072
    if (!task.memory) {
        log.info '[WINDOWMASKER: MK_COUNTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        memory  = (task.memory.toMega()).intValue()
    }

    """
    windowmasker -mk_counts \\
        ${args} \\
        -mem ${memory} \\
        -in ${ref} \\
        -out ${prefix}.txt
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.txt
    """
}
