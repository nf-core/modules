process BISCUIT_QC {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08017a6f4d3c9849d56375068fbe75b8440ed2a1407699958dc3a28759558d1/data':
        'community.wave.seqera.io/library/biscuit:1.8.0.20260217--7d70c1bb73e42ce3' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.txt"), emit: reports
    tuple val("${task.process}"), val('biscuit'), eval("biscuit version |& sed '1!d; s/^.*BISCUIT Version: //'"), emit: versions_biscuit, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_CpGRetentionByReadPos.txt
    touch ${prefix}_CpHRetentionByReadPos.txt
    touch ${prefix}_dup_report.txt
    touch ${prefix}_isize_table.txt
    touch ${prefix}_mapq_table.txt
    touch ${prefix}_strand_table.txt
    touch ${prefix}_totalReadConversionRate.txt
    """
}
