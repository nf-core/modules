process TABIX_TABIX {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), emit: tbi
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    tabix $args $tab

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
