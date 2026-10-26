nextflow.enable.types = true

process RESEQ_ILLUMINAPE {
    tag "${id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/494dd5224f3e7314d36af63a0740e51ecfd5373fe993ab3acf914754870e963a/data' :
        'community.wave.seqera.io/library/reseq:1.1--8536b4889ef6d3b8' }"

    ext [
        prefix: { output_prefix ?: "${id}" },
        args: { reseq_illuminape_args ?: '' }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        bam: Path,
        fasta: Path,
        adapter_fasta: Path?,
        adapter_matrix: Path?,
        reseq_illuminape_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        fastq: files("*.fq.gz").toSorted()
    )

    topic:
    record(process: "${task.process}", tool: 'reseq', version: eval('reseq --version 2>&1 | sed \'s/^.*ReSeq version //\'')) >> 'versions'

    script:
    def adapter_arg = adapter_fasta ? "--adapterFile ${adapter_fasta}" : ''
    def matrix_arg = adapter_matrix ? "--adapterMatrix ${adapter_matrix}" : ''
    """
    reseq \\
        illuminaPE \\
        -j $task.cpus \\
        -r $fasta \\
        -b $bam \\
        $adapter_arg \\
        $matrix_arg \\
        ${task.ext.args} \\
        -1 ${task.ext.prefix}_R1.fq \\
        -2 ${task.ext.prefix}_R2.fq

    gzip -n ${task.ext.prefix}_R1.fq
    gzip -n ${task.ext.prefix}_R2.fq
    """

    stub:
    """
    echo "" | gzip > ${task.ext.prefix}_R1.fq.gz
    echo "" | gzip > ${task.ext.prefix}_R2.fq.gz
    """
}
