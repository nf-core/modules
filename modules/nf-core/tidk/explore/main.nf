process TIDK_EXPLORE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidk:0.2.41--hdbdd923_0':
        'biocontainers/tidk:0.2.41--hdbdd923_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tidk.explore.tsv") , emit: explore_tsv
    tuple val(meta), path("*.top.sequence.txt") , emit: top_sequence, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tidk \\
        explore \\
        $args \\
        $fasta \\
        > ${prefix}.tidk.explore.tsv

    [[ \$(cat ${prefix}.tidk.explore.tsv | wc -l) -gt 1 ]] \\
        && cat \\
        ${prefix}.tidk.explore.tsv \\
        | sed -n 2p \\
        | awk '{print \$1;}' \\
        > ${prefix}.top.sequence.txt \\
        || echo "No sequence identified"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidk: \$(tidk --version | sed 's/tidk //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tidk.explore.tsv
    touch ${prefix}.top.sequence.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidk: \$(tidk --version | sed 's/tidk //')
    END_VERSIONS
    """
}
