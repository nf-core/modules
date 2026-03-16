process PICARD_LIFTOVERVCF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(input_vcf)
    tuple val(meta2), path(dict)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(chain)

    output:
    tuple val(meta), path("*.lifted.vcf.gz"), emit: vcf_lifted
    tuple val(meta), path("*.unlifted.vcf.gz"), emit: vcf_unlifted
    tuple val("${task.process}"), val('picard'), eval("picard LiftoverVcf --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard LiftoverVcf] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        LiftoverVcf \\
        ${args} \\
        --INPUT ${input_vcf} \\
        --OUTPUT ${prefix}.lifted.vcf.gz \\
        --CHAIN ${chain} \\
        --REJECT ${prefix}.unlifted.vcf.gz \\
        --REFERENCE_SEQUENCE ${fasta}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.lifted.vcf.gz
    echo | gzip > ${prefix}.unlifted.vcf.gz
    """
}
