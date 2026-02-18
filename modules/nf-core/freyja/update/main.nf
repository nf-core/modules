process FREYJA_UPDATE {
    tag "$db_name"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:2.0.3--pyhdfd78af_0' :
        'biocontainers/freyja:2.0.3--pyhdfd78af_0' }"

    input:
    val db_name

    output:
    path "${db_name}/*barcodes.*"            , emit: barcodes
    path "${db_name}/*lineages.yml"          , emit: lineages_topology
    path "${db_name}/*pathogen_config.yml"   , emit: config
    path "${db_name}/*curated_lineages.json" , emit: lineages_meta, optional: true
    tuple val("${task.process}"), val('freyja'), eval("freyja --version | sed 's/.* //'"), topic: versions, emit: versions_freyja

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p $db_name
    freyja \\
        update \\
        $args \\
        --outdir $db_name
    """

    stub:
    """
    mkdir $db_name

    touch "${db_name}/usher_barcodes.csv"
    touch "${db_name}/lineages.yml"
    touch "${db_name}/curated_lineages.json"
    touch "${db_name}/pathogen_config.yml"

    """
}
