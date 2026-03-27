process PHARMCAT_MATCHER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b27c134f2226e65c3be9687fdcd6dfb5eebb7998bf1ad89ff396c914fe6d81a/data'
        : 'community.wave.seqera.io/library/pharmcat3:3.1.1--876b7152770ba008'}"

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
