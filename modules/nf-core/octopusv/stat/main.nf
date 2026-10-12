process OCTOPUSV_STAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8378c503810508764badf4504647ea1791fa1f585dc5b77ad5710237609b15ff/data'
:         'community.wave.seqera.io/library/octopusv:1.0.0--5a8c89c65098e801' }"

    input:
    tuple val(meta), path(svcf)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.txt"),  emit: txt,  optional: true
    tuple val(meta), path("*.html"), emit: html, optional: true
    tuple val("${task.process}"), val('octopusv'), eval("python -c \"import importlib.metadata as m; print(m.version('octopusv'))\""), emit: versions_octopusv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def want_json   = args.tokenize().contains('--json')
    def want_report = args.tokenize().contains('--report')

    // Remove output-mode flags because the module handles them below.
    def ext_args = args
        .replaceAll(/(^|\s)--json(?=\s|$)/, ' ')
        .replaceAll(/(^|\s)--report(?=\s|$)/, ' ')
        .trim()

    def json_command = want_json ? """
    octopusv stat \\
        --input-file ${svcf} \\
        --output-file ${prefix}.json \\
        --json \\
        ${ext_args}
    """ : ''

    def text_command = !want_json ? """
    octopusv stat \\
        --input-file ${svcf} \\
        --output-file ${prefix}.txt \\
        ${ext_args} \\
        ${want_report ? '--report' : ''}
    """ : ''

    // JSON and report require separate invocations now
    def json_report_command = want_json && want_report ? """
    octopusv stat \\
        --input-file ${svcf} \\
        --output-file ${prefix}.txt \\
        ${ext_args} \\
        --report

    rm ${prefix}.txt
    """ : ''

    """
    ${json_command}
    ${text_command}
    ${json_report_command}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def want_json   = args.tokenize().contains('--json')
    def want_report = args.tokenize().contains('--report')

    def json_output = want_json ? "touch ${prefix}.json" : ''
    def text_output = !want_json ? "touch ${prefix}.txt" : ''
    def html_output = want_report ? "touch ${prefix}.html" : ''

    """
    echo ${args}

    ${json_output}
    ${text_output}
    ${html_output}
    """
}
