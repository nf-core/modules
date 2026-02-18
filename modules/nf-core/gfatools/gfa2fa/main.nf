process GFATOOLS_GFA2FA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfatools:0.5--h577a1d6_5':
        'biocontainers/gfatools:0.5--h577a1d6_5' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val("${task.process}"), val('gfatools'), eval("gfatools version | sed '1!d; s/.* //'"), topic: versions, emit: versions_gfatools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gfatools \\
        gfa2fa \\
        $args \\
        $gfa \\
        | gzip -c > ${prefix}.fasta.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.fasta.gz
    """
}
