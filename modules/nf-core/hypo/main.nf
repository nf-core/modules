// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HYPO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13780867e01535dbfb55bc66e19de88b396d9c21d1077ef52f7134f7c14253c3/data'
        : 'community.wave.seqera.io/library/hypo:1.0.3--6e4b00f92bd33cc2'}"

    input:
    tuple val(meta), path(sr_bam)
    tuple val(meta2), path(reads)
    path draft
    val genome_size
    val reads_coverage

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val("hypo"), val("1.0.3"), topic: versions, emit: versions_hypo

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${reads.join(' ')}" | tr " " "\\n" > sr.fofn
    hypo \\
        -r @sr.fofn \\
        -d ${draft} \\
        -b ${sr_bam} \\
        -c ${reads_coverage} \\
        -s ${genome_size} \\
        -t ${task.cpus} \\
        -o ${prefix}.fasta \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    """
}
