process GATK4_UNMARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/927ff9bb80d65b425cbe752db6648a84043feff6e8ca90e60f9ff6ddbe8938d5/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:c1e4292d6ee27439'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_UnmarkDuplicates"
    def input_list = bam.collect { "--input ${it}" }.join(' ')
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK UnmarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        UnmarkDuplicates \\
        ${input_list} \\
        --output ${prefix}.bam \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_UnmarkDuplicates"
    prefix_no_suffix = task.ext.prefix ? prefix.tokenize('.')[0] : "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
