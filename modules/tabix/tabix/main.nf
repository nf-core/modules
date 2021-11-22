process TABIX_TABIX {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"
    }

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), emit: tbi
    path  "versions.yml"          , emit: versions

    script:
    """
    tabix $options.args $tab

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
