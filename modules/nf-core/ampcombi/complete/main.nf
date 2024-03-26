process AMPCOMBI_COMPLETE {
    tag "ampcombi"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.2.2--pyhdfd78af_0':
        'biocontainers/ampcombi:0.2.2--pyhdfd78af_0' }"

    input:
    path(summaries)

    output:
    tuple val(meta), path("Ampcombi_summary.tsv")   , emit: tsv
    tuple val(meta), path("Ampcombi_complete.log")  , optional:true, emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def summary_input = summaries.isDirectory() ? "--summaries_directory ${summaries}/" : "--summaries_files '${summaries.collect{"$it"}.join("' '")}'"
    """
    ampcombi complete \\
        $args \\
        ${summary_input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def summary_input = summaries.isDirectory() ? "--summaries_directory ${summaries}/" : "--summaries_files '${summaries.collect{"$it"}.join("' '")}'"
    """
    touch Ampcombi_summary.tsv
    touch Ampcombi_complete.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: \$(ampcombi --version | sed 's/ampcombi //')
    END_VERSIONS
    """
}
