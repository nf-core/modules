process PHARMCAT_MATCHER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7dd711a2b130b55d33e119a346ef8040191bf7834a3c393ed6e29d7d9026d5e/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.2.0--5126bb296d1e59ac'}"

    input:
    tuple val(meta), path(vcf), path(index)
    val genes

    output:
    tuple val(meta), path("*.match.json")                                                                       ,   emit: matcher_json
    tuple val(meta), path("*.match.html")                                                     , optional: true  ,   emit: matcher_html
    tuple val("${task.process}"), val('pharmcat'), eval("pharmcat --version | cut -f2 -d ' '"), topic: versions ,   emit: versions_pharmcat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.pharmcat"

    // Genes lists
    def genes_join = genes instanceof List ? genes.collect().join(',') : null
    def genes_cmd  = genes_join ? "--genes ${genes_join}" : ""

    """
    pharmcat \\
        -vcf ${vcf} \\
        --base-filename ${prefix} \\
        --output-dir . \\
        -matcher \\
        --samples ${meta.id} \\
        ${args} \\
        ${genes_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.pharmcat"
    """
    echo ${args} >/dev/null

    touch ${prefix}.match.json
    touch ${prefix}.match.html
    """
}
