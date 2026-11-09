// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMTOOLS_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8345056363cc2598d6185cfab19d67eda17a5ea73a4f9b56c25322347b97aef9/data'
        : 'community.wave.seqera.io/library/htslib_samtools_xz:ca51bc8a00d53523'}"

    input:
    tuple val(meta), path(input)
    val interleave

    output:
    tuple val(meta), path("*_{1,2}.fasta.gz"), emit: fasta, optional: true
    tuple val(meta), path("*_interleaved.fasta.gz"), emit: interleaved, optional: true
    tuple val(meta), path("*_singleton.fasta.gz"), emit: singleton, optional: true
    tuple val(meta), path("*_other.fasta.gz"), emit: other, optional: true
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = interleave && !meta.single_end
        ? "| bgzip > ${prefix}_interleaved.fasta.gz"
        : meta.single_end
            ? "-1 ${prefix}_1.fasta.gz -s ${prefix}_singleton.fasta.gz"
            : "-1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz -s ${prefix}_singleton.fasta.gz"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        fasta \\
        ${args} \\
        --threads ${task.cpus - 1} \\
        -0 ${prefix}_other.fasta.gz \\
        ${input} \\
        ${output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputs = []
    if (interleave && !meta.single_end) {
        outputs << "echo | bgzip -c > ${prefix}_interleaved.fasta.gz"
    }
    else if (meta.single_end) {
        outputs << "echo | bgzip -c > ${prefix}_1.fasta.gz"
        outputs << "echo | bgzip -c > ${prefix}_singleton.fasta.gz"
    }
    else {
        outputs << "echo | bgzip -c > ${prefix}_1.fasta.gz"
        outputs << "echo | bgzip -c > ${prefix}_2.fasta.gz"
        outputs << "echo | bgzip -c > ${prefix}_singleton.fasta.gz"
    }
    outputs << "echo | bgzip -c > ${prefix}_other.fasta.gz"

    """
    ${outputs.join('\n')}
    """
}
