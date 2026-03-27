process PICARD_MERGESAMFILES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data'
        : 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('picard'), eval("picard MergeSamFiles --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_files = bams.sort()
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    if (bam_files.size() > 1) {
        """
        picard \\
            -Xmx${avail_mem}M \\
            MergeSamFiles \\
            ${args} \\
            ${'--INPUT ' + bam_files.join(' --INPUT ')} \\
            --OUTPUT ${prefix}.bam
        """
    }
    else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
