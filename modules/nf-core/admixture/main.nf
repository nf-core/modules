
process ADMIXTURE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::admixture=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/admixture:1.3.0--0':
        'quay.io/biocontainers/admixture:1.3.0--0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val K


    output:
    tuple val(meta), path("*.Q")    , emit: ancestry_fractions
    tuple val(meta), path("*.P")    , emit: allele_frequencies
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    admixture \\
        $bed \\
        $K \\
        -j$task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(echo \$(admixture 2>&1) | head -n 1  )
    END_VERSIONS

    """
}
