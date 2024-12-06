process AMPCOMBI2_COMPLETE {
    tag "ampcombi2"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
    path(summaries)

    output:
    path("Ampcombi_summary.tsv") , emit: tsv
    path("Ampcombi_complete.log"), emit: log, optional:true
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ampcombi complete \\
        --summaries_files '${summaries.collect{"$it"}.join("' '")}' \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch Ampcombi_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
