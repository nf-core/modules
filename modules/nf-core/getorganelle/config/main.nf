process GETORGANELLE_CONFIG {
    tag "${organelle_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0'
        : 'quay.io/biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0'}"

    input:
    val(organelle_type)

    output:
    tuple val(organelle_type), path("getorganelle"), emit: db
    tuple val("${task.process}"), val('getorganelle'), eval("get_organelle_config.py --version |& sed 's/^GetOrganelle v//'"), topic: versions, emit: versions_getorganelle

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    get_organelle_config.py \\
        ${args} \\
        -a ${organelle_type} \\
        --config-dir getorganelle
    """

    stub:
    """
    mkdir -p getorganelle/{LabelDatabase,SeedDatabase}
    touch getorganelle/{LabelDatabase,SeedDatabase}/${organelle_type}.fasta
    """
}
