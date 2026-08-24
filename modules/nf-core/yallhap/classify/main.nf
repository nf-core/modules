process YALLHAP_CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yallhap:0.4.0--pyhdfd78af_0':
        'biocontainers/yallhap:0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(tree)
    tuple val(meta3), path(snp_db)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.tsv") , emit: tsv , optional: true
    tuple val("${task.process}"), val('yallhap'), eval("yallhap --version | sed 's/.*version //'"), emit: versions_yallhap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = task.ext.format ?: 'json'
    def output_ext = format == 'tsv' ? 'tsv' : 'json'
    """
    yallhap \\
        classify \\
        ${vcf} \\
        --tree ${tree} \\
        --snp-db ${snp_db} \\
        --format ${format} \\
        --output ${prefix}.${output_ext} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = task.ext.format ?: 'json'
    def output_ext = format == 'tsv' ? 'tsv' : 'json'
    """
    touch ${prefix}.${output_ext}
    """
}
