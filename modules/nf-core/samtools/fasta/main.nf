process SAMTOOLS_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

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
        ? "| gzip > ${prefix}_interleaved.fasta.gz"
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
        outputs << "echo | gzip > ${prefix}_interleaved.fasta.gz"
    }
    else if (meta.single_end) {
        outputs << "echo | gzip > ${prefix}_1.fasta.gz"
        outputs << "echo | gzip > ${prefix}_singleton.fasta.gz"
    }
    else {
        outputs << "echo | gzip > ${prefix}_1.fasta.gz"
        outputs << "echo | gzip > ${prefix}_2.fasta.gz"
        outputs << "echo | gzip > ${prefix}_singleton.fasta.gz"
    }
    outputs << "echo | gzip > ${prefix}_other.fasta.gz"

    """
    ${outputs.join('\n')}
    """
}
