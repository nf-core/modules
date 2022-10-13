process ENDORSPY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::endorspy=0.4" : null)
        'https://depot.galaxyproject.org/singularity/endorspy:0.4--hdfd78af_0':
        "${params.docker_registry ?: 'quay.io/biocontainers'}/endorspy:0.4--hdfd78af_0" }

    input:
    tuple val(meta), path(stats), path(stats_optional)

    output:
    tuple val(meta), path("*_mqc.json"), emit: json
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def optionalstats = stats_optional ? "statsOptional=${stats_optional}" : ''

    """
    endorspy \\
        $stats \\
        $optionalstats \\
        $args \\
        -o json \\
        -n $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        endorspy: \$(echo \$(endorspy --version 2>&1) | sed 's/^endorS.py //' )
    END_VERSIONS
    """
}
