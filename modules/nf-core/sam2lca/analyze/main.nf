
process SAM2LCA_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.4--pyhdfd78af_0':
        'quay.io/biocontainers/sam2lca:1.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(database)

    output:
    tuple val(meta), path("*.csv")  , emit: csv
    tuple val(meta), path("*.json") , emit: json
    tuple val(meta), path("*.bam")  , emit: bam     , optional: true
    tuple val("${task.process}"), val('sam2lca'), eval("sam2lca --version | sed 's/.*version //'"), topic: versions, emit: versions_sam2lca

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    touch ${prefix}.json
    """
}
