process BISMARK_SUMMARY {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bismark=0.24.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.0--hdfd78af_0' :
        'quay.io/biocontainers/bismark:0.24.0--hdfd78af_0' }"

    input:
    val(bam)
    path(align_report)
    path(dedup_report)
    path(splitting_report)
    path(mbias)

    output:
    path  "*report.{html,txt}", emit: summary
    path  "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark2summary ${bam.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
