process PICARD_CLEANSAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.27.4" : null)
    def container_image = "/picard:2.27.4--hdfd78af_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard CleanSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CleanSam  \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CleanSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
