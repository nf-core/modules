process DIAMOND_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63594b00a9e68e690f398d9b9788a883753442dedcce590b4ba5aa97cff200f1/data'
        : 'community.wave.seqera.io/library/diamond:2.2.2--910ae987965d0d4d'}"

    input:
    tuple val(meta), path(db)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('diamond'), eval("diamond --version | sed 's/diamond version //g'"), emit: versions_diamond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.memory.toKilo() + 'K'
    def memarg = "-M ${mem}"
    """
    diamond \\
        cluster \\
        ${args} \\
        ${memarg} \\
        -p ${task.cpus} \\
        -d ${db} \\
        -o ${prefix}.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    touch ${prefix}.tsv
    """
}
