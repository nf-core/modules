process VSEARCH_SINTAX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/eec4fdf832cbd591111a702a1942476674a71e76fafbd5ae3fc18018b9ebe2e2/data'
:         'community.wave.seqera.io/library/vsearch:2.31.0--80e3b2da3ae2f3d8' }"

    input:
    tuple val(meta), path(queryfasta)
    path db

    output:
    tuple val(meta), path('*.tsv')   , optional: true, emit: tsv
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | sed -n "1s/.*v\\([0-9.]*\\).*/\\\\1/p"'), emit: versions_vsearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vsearch \\
        --sintax $queryfasta \\
        --db $db \\
        --threads $task.cpus \\
        $args \\
        --tabbedout ${prefix}.tsv
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
