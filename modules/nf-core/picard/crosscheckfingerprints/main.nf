process PICARD_CROSSCHECKFINGERPRINTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(input1), path(input1_index), path(input2), path(input2_index), path(haplotype_map)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.crosscheck_metrics.txt"), emit: crosscheck_metrics
    tuple val("${task.process}"), val('picard'), eval("picard CrosscheckFingerprints --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input1_cmd = input1.collect { f1 -> "--INPUT ${f1}" }.join(' ')
    def input2_cmd = input2.collect { f2 -> "--SECOND_INPUT ${f2}" }.join(' ')
    def reference_cmd = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard CrosscheckFingerprints] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CrosscheckFingerprints \\
        ${input1_cmd} \\
        ${input2_cmd} \\
        ${reference_cmd} \\
        --HAPLOTYPE_MAP ${haplotype_map} \\
        --OUTPUT ${prefix}.crosscheck_metrics.txt \\
        --NUM_THREADS ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.crosscheck_metrics.txt
    """
}
