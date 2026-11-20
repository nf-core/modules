nextflow.enable.types = true

process RESEQ_ILLUMINAPE {
    tag "${id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/49/494dd5224f3e7314d36af63a0740e51ecfd5373fe993ab3acf914754870e963a/data' :
        'community.wave.seqera.io/library/reseq:1.1--8536b4889ef6d3b8' }"

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
        fastq1: file("*_R1.fastq.gz"),
        fastq2: file("*_R2.fastq.gz")
    )

    topic:
    record(process: "${task.process}", tool: 'reseq', version: eval('reseq --version 2>&1 | sed \'s/^.*ReSeq version //\'')) >> 'versions'

    script:
    def ext = task.ext
    def args = reseq_illuminape_args ?: ext.args ?: ''
    def prefix = output_prefix ?: ext.prefix ?: id
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
        $args \\
        -1 ${prefix}_R1.fastq \\
        -2 ${prefix}_R2.fastq

    gzip -n ${prefix}_R1.fastq
    gzip -n ${prefix}_R2.fastq
    """

    stub:
    def ext = task.ext
    def prefix = output_prefix ?: ext.prefix ?: id
    """
    echo "" | gzip > ${prefix}_R1.fastq.gz
    echo "" | gzip > ${prefix}_R2.fastq.gz
    """
}
