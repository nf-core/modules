process SEQKIT_SANA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}${extension}"), emit: reads
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = reads.getName() - reads.getSimpleName()
    """
    seqkit sana \\
        ${args} \\
        ${reads} \\
        -o ${prefix}${extension} > ${prefix}.log 2>&1
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = reads.getName() - reads.getSimpleName()
    def create_cmd = extension.endsWith('.gz') ? "echo '' | gzip >" : "touch"
    """
    echo ${args}

    ${create_cmd} ${prefix}${extension}
    touch ${prefix}.log
    """
}
