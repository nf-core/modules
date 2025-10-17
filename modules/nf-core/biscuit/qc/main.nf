process BISCUIT_QC {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33a9ca30b4154f11253c8d91a75382065dcb8282ba99b74dbee59ed8faceabd7/data':
        'community.wave.seqera.io/library/biscuit:1.5.0.20240506--ca92d9d0a37b5fa8' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.txt"), emit: reports
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def se = meta.single_end ? "-s" : ""
    """
    ln -sf \$(readlink $fasta) $index/$fasta

    biscuit qc \\
        $args \\
        $se \\
        $index/$fasta \\
        $bam \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_CpGRetentionByReadPos.txt
    touch ${prefix}_CpHRetentionByReadPos.txt
    touch ${prefix}_dup_report.txt
    touch ${prefix}_isize_table.txt
    touch ${prefix}_mapq_table.txt
    touch ${prefix}_strand_table.txt
    touch ${prefix}_totalReadConversionRate.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}
