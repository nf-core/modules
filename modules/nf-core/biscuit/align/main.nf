process BISCUIT_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/5021d3f67d0221a5fa1eb692573a20fa28a5a552754d89c4351eb0eba42eadb0/data':
        'community.wave.seqera.io/library/biscuit_samtools:1e5147589c9d9ec1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    tuple val("${task.process}"), val('biscuit'), eval("biscuit version |& sed '1!d; s/^.*BISCUIT Version: //'"), emit: versions_biscuit, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version |& sed '1!d; s/^.*samtools //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = [(task.cpus * 0.9) as int, 1].max()
    def samtools_cpus = (task.cpus - biscuit_cpus < 1) ? biscuit_cpus : (task.cpus - biscuit_cpus)
    """
    ln -sf \$(readlink $fasta) $index/$fasta

    biscuit align \\
        $args \\
        -@ $biscuit_cpus \\
        $index/$fasta \\
        $reads \\
        | samtools sort $args2 --threads $samtools_cpus --write-index -o ${prefix}.bam -O BAM -

    samtools index ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}
