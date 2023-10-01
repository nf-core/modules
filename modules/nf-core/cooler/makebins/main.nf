process COOLER_MAKEBINS {
    tag "${meta.id}}"
    label 'process_low'

    conda "bioconda::cooler=0.9.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.9.2--pyh7cba7a3_0' :
        'biocontainers/cooler:0.9.2--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(chromsizes), val(cool_bin)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooler makebins \\
        $args \\
        ${chromsizes} \\
        ${cool_bin} > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
