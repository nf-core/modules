process NEXTCLADE_DATASETGET {
    tag "$dataset"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::nextclade=1.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:1.9.0--h9ee0642_0' :
        'quay.io/biocontainers/nextclade:1.9.0--h9ee0642_0' }"

    input:
    val dataset

    output:
    path "$prefix"     , emit: dataset
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dataset}"
    """
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name $dataset \\
        --output-dir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1)
    END_VERSIONS
    """
}
