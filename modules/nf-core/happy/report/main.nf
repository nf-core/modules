process HAPPY_REPORT {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6e1aa20442ec807e9aa59cc5b46811e6758667c295610a745c5ab9cdb23c4e0d/data':
        'community.wave.seqera.io/library/pip_ga4gh-happy-report:72e6bf058d448a4c' }"

    input:
    tuple val(meta), val(csv_meta), path(csv, arity: '1..*')

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val("${task.process}"), val('ga4gh-happy-report'), val('0.1.2'), topic: versions, emit: versions_happy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def csv_list = csv instanceof List ? csv : [csv]
    def csv_meta_list = csv_meta instanceof List ? csv_meta : [csv_meta]
    assert csv_meta_list.size() == csv_list.size() : "HAPPY_REPORT: 'csv_meta' and 'csv' must be the same length (got ${csv_meta_list.size()} csv_meta entries, ${csv_list.size()} files)"
    def labelled_csv = [csv_meta_list, csv_list]
        .transpose()
        .collect { m, file -> m.comparison_method ? "${m.method}_${m.comparison_method}:${file}" : "${m.method}:${file}" }
        .join(' ')
    """
    rep.py \\
        $args \\
        $labelled_csv \\
        -o ${prefix}.html
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    """
}
