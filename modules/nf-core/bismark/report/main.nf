process BISMARK_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    def container_image = "bismark:0.23.0--0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


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
    bismark2report $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
