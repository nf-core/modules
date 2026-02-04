
process SAM2LCA_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.4--pyhdfd78af_0':
        'biocontainers/sam2lca:1.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(database)

    output:
    tuple val(meta), path("*.csv")  , emit: csv
    tuple val(meta), path("*.json") , emit: json
    tuple val(meta), path("*.bam")  , emit: bam     , optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def make_db = database ? "" : "mkdir sam2lca_db"
    def database_path = database ? "${database}" : "sam2lca_db"
    """
    $make_db
    sam2lca \\
        -d $database_path \\
        analyze \\
        $args \\
        -o ${prefix} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(sam2lca --version | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(sam2lca --version | sed 's/.*version //')
    END_VERSIONS
    """
}
