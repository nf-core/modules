process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::svdb=2.7.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.7.1--py27h4329609_0':
        'quay.io/biocontainers/svdb:2.7.1--py27h4329609_0' }"

    input:
    tuple val(meta), path(vcf)
    path (vcf_dbs)

    output:
    tuple val(meta), path("*_query.vcf"), emit: vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    svdb \\
        --query \\
        --query_vcf $vcf \\
        --db ${vcf_dbs.join(',')} \\
        --prefix ${prefix}
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_query.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """
}
