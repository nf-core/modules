process ENDORSPY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::endorspy=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/endorspy:1.3--hdfd78af_0':
        'biocontainers/endorspy:1.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(stats_raw), path(stats_qualityfiltered), path(stats_deduplicated)

    output:
    tuple val(meta), path("*_mqc.json"), emit: json
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def optionalraw = stats_raw ? "-r ${stats_raw}" : ''
    def optionalqualityfiltered = stats_qualityfiltered ? "-q ${stats_qualityfiltered}" : ''
    def optionaldeduplicated = stats_deduplicated ? "-d ${stats_deduplicated}" : ''

    if ( stats_qualityfiltered && !stats_raw && !stats_deduplicated ) error "ERROR: only input channel stats_qualityfiltered provided. No stats can be calculated. Add at least one additional input channel: stats_raw or stats_deduplicated"


    """
    endorspy \\
        $optionalraw \\
        $optionalqualityfiltered \\
        $optionaldeduplicated \\
        $args \\
        -o json \\
        -n $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        endorspy: \$(echo \$(endorspy --version 2>&1) | sed 's/^endorS.py //' )
    END_VERSIONS
    """
}
