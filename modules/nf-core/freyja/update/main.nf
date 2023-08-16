process FREYJA_UPDATE {
    tag "$db_name"
    label 'process_single'

    conda "bioconda::freyja=1.3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freyja:1.3.12--pyhdfd78af_0':
        'biocontainers/freyja:1.3.12--pyhdfd78af_0' }"

    input:
    val db_name

    output:
    path "${db_name}/usher_barcodes.csv"   , emit: barcodes
    path "${db_name}/lineages.yml"         , emit: lineages_topology
    path "${db_name}/curated_lineages.json", emit: lineages_meta
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p $db_name
    freyja \\
        update \\
        --outdir $db_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """

    stub:
    """
    mkdir $db_name

    touch "${db_name}/usher_barcodes.csv"
    touch "${db_name}/lineages.yml"
    touch "${db_name}/curated_lineages.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
    END_VERSIONS
    """
}
