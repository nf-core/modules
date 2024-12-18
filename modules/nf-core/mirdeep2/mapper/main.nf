process MIRDEEP2_MAPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mirdeep2:2.0.1.2--0':
        'biocontainers/mirdeep2:2.0.1.2--0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index, stageAs: '*')

    output:
    tuple val(meta), path("*.fa"), path("*.arf"), emit: outputs
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.0.1'

    """
    mapper.pl \\
    ${reads} \\
    $args \\
    -p ${index}/${meta2.id}  \\
    -s ${prefix}_collapsed.fa \\
    -t ${prefix}_reads_collapsed_vs_${meta2.id}_genome.arf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirdeep2: \$(echo "$VERSION")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.0.1'
    """
    touch ${prefix}.fa
    touch ${prefix}reads_vs_refdb.arf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirdeep2: \$(echo "$VERSION")
    END_VERSIONS
    """
}
