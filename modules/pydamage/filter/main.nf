process PYDAMAGE_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pydamage=0.70" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pydamage:0.70--pyhdfd78af_0' :
        'quay.io/biocontainers/pydamage:0.70--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("pydamage_results/pydamage_filtered_results.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    pydamage \\
        filter \\
        $args \\
        $csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pydamage: \$(echo \$(pydamage --version 2>&1) | sed -e 's/pydamage, version //g')
    END_VERSIONS
    """
}
