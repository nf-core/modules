process BISMARK_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.2--hdfd78af_0' :
        'biocontainers/bismark:0.24.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(align_report), path(dedup_report), path(splitting_report), path(mbias)

    output:
    tuple val(meta), path("*report.{html,txt}"), emit: report
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark2report ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.txt
    touch ${prefix}.report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
