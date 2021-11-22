process PYDAMAGE_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pydamage=0.62" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pydamage:0.62--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/pydamage:0.62--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("pydamage_results/pydamage_filtered_results.csv"), emit: csv
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """

    pydamage \\
        filter \\
        $args \\
        $csv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(pydamage --version 2>&1) | sed -e 's/pydamage, version //g')
    END_VERSIONS
    """
}
