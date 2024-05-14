
process SAM2LCA_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sam2lca=1.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.2--pyhdfd78af_1':
        'biocontainers/sam2lca:1.1.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(database)

    output:
    tuple val(meta), path("*.csv") ,                 emit: csv
    tuple val(meta), path("*.json"),                 emit: json
    tuple val(meta), path("*.bam") , optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def make_db = database ? "" : "mkdir sam2lca_db"
    def database = database ? "${database}" : "sam2lca_db"
    """
    $make_db
    sam2lca \\
        -d $database \\
        analyze \\
        $args \\
        -o ${prefix} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/^sam2lca, version //' ))
    END_VERSIONS
    """
}
