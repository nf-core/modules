process PHARMCAT_REPORTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.1.1--876b7152770ba008'}"

    input:
    tuple val(meta), path(phenotypes)

    output:
    tuple val(meta), path("*.report.json")                                                    , optional: true  ,   emit: report_json
    tuple val(meta), path("*.report.html")                                                    , optional: true  ,   emit: report_html
    tuple val(meta), path("*.report.tsv")                                                     , optional: true  ,   emit: report_tsv
    tuple val("${task.process}"), val('pharmcat'), eval("pharmcat --version | cut -f2 -d ' '"), topic: versions ,   emit: versions_pharmcat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}.pharmcat"

    """
    pharmcat \\
        -reporter \\
        --reporter-input ${phenotypes} \\
        --base-filename ${prefix} \\
        --output-dir . \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.pharmcat"
    """
    echo ${args} >/dev/null

    touch ${prefix}.report.json
    touch ${prefix}.report.html
    touch ${prefix}.report.tsv
    """
}
