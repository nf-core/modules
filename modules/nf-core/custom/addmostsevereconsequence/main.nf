process CUSTOM_ADDMOSTSEVERECONSEQUENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbeb20c898a76bec809629320ece9e1f84a3e355e96568bfbe14b9f411bdf3e7/data':
        'community.wave.seqera.io/library/htslib_python:9c6265e98ef06930' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(variant_consequences)


    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('addmostsevereconsequence'), val("1.1"), topic: versions, emit: versions_addmostsevereconsequence
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    add_most_severe_consequence.py \\
        $args \\
        --file_in $vcf \\
        --file_out ${prefix}.vcf \\
        --variant_csq $variant_consequences

    bgzip \\
        $args2 \\
        --threads ${task.cpus} \\
        ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    """
}
