process PHARMCAT_REPORTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7dd711a2b130b55d33e119a346ef8040191bf7834a3c393ed6e29d7d9026d5e/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.2.0--5126bb296d1e59ac'}"

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
