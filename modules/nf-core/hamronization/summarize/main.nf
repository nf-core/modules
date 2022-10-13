process HAMRONIZATION_SUMMARIZE {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::hamronization=1.1.1" : null)
    def container_image = "/hamronization:1.1.1--pyhdfd78af_0"
                                                          container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path(reports)
    val(format)

    output:
    path("hamronization_combined_report.json"), optional: true, emit: json
    path("hamronization_combined_report.tsv") , optional: true, emit: tsv
    path("hamronization_combined_report.html"), optional: true, emit: html
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def outformat = format == 'interactive' ? 'html' : format
    """
    hamronize \\
        summarize \\
        ${reports.join(' ')} \\
        -t ${format} \\
        $args \\
        -o hamronization_combined_report.${outformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(echo \$(hamronize --version 2>&1) | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
