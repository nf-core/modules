// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CNVKIT_PURITY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/3935f8d7507f85fd20e073f0625bb6bdff8aa1b6044f4bbb720c9a8063ea0390/data'
        : 'community.wave.seqera.io/library/cnvkit:0.9.14--288e98d6210b7304'}"

    input:
    tuple val(meta), path(cns), path(vcf)

    output:
    tuple val(meta), path("*.purity.tsv"), emit: purity
    tuple val("${task.process}"), val('cnvkit'), eval("cnvkit.py version | sed -e 's/cnvkit v//g'"), topic: versions, emit: versions_cnvkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_cmd = vcf ? "-v ${vcf}" : ""

    """
    cnvkit.py purity \\
        ${cns} \\
        ${args} \\
        ${vcf_cmd} \\
        --output ${prefix}.purity.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.purity.tsv
    """
}
