process HAMRONIZATION_SUMMARIZE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hamronization:1.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/hamronization:1.1.9--pyhdfd78af_0'}"

    input:
    path reports
    val format

    output:
    path ("hamronization_combined_report.json"), optional: true, emit: json
    path ("hamronization_combined_report.tsv"), optional: true, emit: tsv
    path ("hamronization_combined_report.html"), optional: true, emit: html
    tuple val("${task.process}"), val('hamronization'), eval("hamronize --version 2>&1 | sed 's/hamronize //'"), topic: versions, emit: versions_hamronization

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
        ${args} \\
        -o hamronization_combined_report.${outformat}
    """

    stub:
    def outformat = format == 'interactive' ? 'html' : format
    """
    touch hamronization_combined_report.${outformat}
    """
}
