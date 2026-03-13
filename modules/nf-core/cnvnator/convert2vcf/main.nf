process CNVNATOR_CONVERT2VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a54de0a286dcc6f73972be37104995ba96e37ef404f9d2bff294a7e8be5b7a92/data':
        'community.wave.seqera.io/library/cnvnator:0.4.1--5a467cfadbbc668d' }"

    input:
    tuple val(meta), path(calls)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('cnvnator'), eval("cnvnator 2>&1 | sed -n '3s/CNVnator v//p'"), topic: versions, emit: versions_cnvnator

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    cnvnator2VCF.pl \\
        ${calls} \\
        $args \\
        > ${prefix}.vcf
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.vcf
    """
}
