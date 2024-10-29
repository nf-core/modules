process PICARD_POSITIONBASEDDOWNSAMPLESAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0' :
        'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), val(fraction)

    output:
    tuple val(meta), path("*.ds*.bam")        , emit: bam
    tuple val(meta), path("*.ds*.bai")        , emit: bai, optional:true
    tuple val(meta), env(ACTUAL_NUM_READS)    , emit: num_reads, optional:true
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard PositionBasedDownsampleSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        PositionBasedDownsampleSam \\
        $args \\
        --CREATE_INDEX \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.ds.bam \\
        --FRACTION ${fraction} 2> tool_stderr

    ACTUAL_NUM_READS=\$(tail -n 10 tool_stderr | grep Kept | sed -E 's/.*Kept ([0-9]+) out of.*/\\1/')
    mv "${prefix}.ds.bam" "${prefix}.ds\${ACTUAL_NUM_READS}.bam"
    mv "${prefix}.ds.bai" "${prefix}.ds\${ACTUAL_NUM_READS}.bai"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard PositionBasedDownsampleSam --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ds10.bam
    touch ${prefix}.ds10.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard PositionBasedDownsampleSam --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
