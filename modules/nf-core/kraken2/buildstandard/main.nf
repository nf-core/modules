process KRAKEN2_BUILDSTANDARD {
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7d3a961582e778e2aa43a61d61a09c7fe1871808f3d9413dcf22102bdecf09c/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:d1f7e5e385456c6c' }"

    input:
    val cleaning

    output:
    path("$prefix"),     emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "kraken2_standard_db"
    runclean = cleaning ? "kraken2-build --clean --db ${db}" : ""
    """
    kraken2-build \\
        --standard \\
        $args \\
        --threads ${task.cpus} \\
        --db ${prefix}
    $runclean
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "kraken2_standard_db"
    """
    mkdir -p "$prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
