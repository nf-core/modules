process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::svdb=2.8.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.8.1--py39h5371cbf_0':
        'quay.io/biocontainers/svdb:2.8.1--py39h5371cbf_0' }"

    input:
    tuple val(meta), path(vcf)
    path (vcf_dbs)

    output:
    tuple val(meta), path("*_query.vcf")    , emit: vcf
    path "versions.yml"                     , emit: versions

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
