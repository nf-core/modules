process BISMARK_COVERAGE2CYTOSINE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bismark=0.23.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.0--hdfd78af_0' :
        'biocontainers/bismark:0.24.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(coverage_file)
    path index

    output:
    tuple val(meta), path("*.cov.gz"), optional: true      , emit: coverage
    tuple val(meta), path("*report.txt.gz")                , emit: report
    tuple val(meta), path("*cytosine_context_summary.txt") , emit: summary
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    coverage2cytosine \\
        $coverage_file \\
        --genome $index \\
        --output ${prefix} \\
        --gzip \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
