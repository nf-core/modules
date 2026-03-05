process RESEQ_ILLUMINAPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/reseq:1.1--ff843e111e12fa4e' :
        'community.wave.seqera.io/library/reseq:1.1--8536b4889ef6d3b8' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.fq.gz"), emit: fastq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reseq \\
        illuminaPE \\
        -j $task.cpus \\
        -r $fasta \\
        -b $bam \\
        $args \\
        -1 ${prefix}_R1.fq \\
        -2 ${prefix}_R2.fq

    gzip ${prefix}_R1.fq
    gzip ${prefix}_R2.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reseq: \$(reseq --version 2>&1 | sed 's/^.*ReSeq version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}_R1.fq.gz
    echo '' | gzip > ${prefix}_R2.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reseq: \$(reseq --version 2>&1 | sed 's/^.*ReSeq version: //')
    END_VERSIONS
    """
}
