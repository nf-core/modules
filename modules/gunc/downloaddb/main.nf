process GUNC_DOWNLOADDB {
    tag '$db_name'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gunc=1.0.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gunc:1.0.5--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gunc:1.0.5--pyhdfd78af_0"
    }

    input:
    val db_name

    output:
    path "*.dmnd"       , emit: db
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    gunc download_db . -db $db_name $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( gunc --version )
    END_VERSIONS
    """
}
