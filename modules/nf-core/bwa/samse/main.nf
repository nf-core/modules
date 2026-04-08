process BWA_SAMSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(reads), path(sai)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), topic: versions, emit: versions_bwa
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ? "-r ${meta.read_group}" : ""

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa samse \\
        $args \\
        $read_group \\
        \$INDEX \\
        $sai \\
        $reads | samtools sort -@ ${task.cpus - 1} -O bam - > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
