process ABRA2 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/abra2:2.24--hdcf5f25_3'
        : 'biocontainers/abra2:2.24--hdcf5f25_3'}"

    input:
    tuple val(meta), path(bams), path(bais)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(targets)
    tuple val(meta5), path(gtf)
    tuple val(meta6), path(known_indels)

    output:
    tuple val(meta), path("*.bam"),     emit: bam
    tuple val(meta), path("*.bam.bai"), emit: bai, optional: true
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def input_bams       = bams.join(",")
    def targets_arg      = targets         ? "--targets ${targets}"     : ""
    def gtf_arg          = gtf             ? "--gtf ${gtf}"             : ""
    def known_indels_arg = known_indels    ? "--in-vcf ${known_indels}" : ""
    def output_str       = bams.collect { bam -> "${prefix}.${bam.name}" }.join(",")
    def memory_cmd       = task.memory     ? "export JAVA_TOOL_OPTIONS='-Xmx${task.memory.toGiga()}G'" : ""
    """
    ${memory_cmd}
    abra2 \\
        --in ${input_bams} \\
        --out ${output_str} \\
        --ref ${fasta} \\
        --threads ${task.cpus} \\
        ${targets_arg} \\
        ${gtf_arg} \\
        ${known_indels_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: \$(abra2 2>&1 | grep 'Abra version:' | sed 's/.*Abra version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.abra.bam
    touch ${prefix}.abra.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: \$(abra2 2>&1 | grep 'Abra version:' | sed 's/.*Abra version: //')
    END_VERSIONS
    """
}
