process PBMM2_ALIGN {
    tag "$meta.id"
    label 'process_high'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/86/8678f70f464a6afd39f737c492a6e7d627ee7bb8a9cf54fe5008b564834084b5/data':
        'community.wave.seqera.io/library/pbmm2:26.2.0--37598eea709c00f6' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('pbmm2'), eval("pbmm2 --version  | sed '1s/pbmm2 //;q "), topic: versions, emit: versions_pbmm2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmm2 \\
        align \\
        $args \\
        $fasta \\
        $bam \\
        ${prefix}.bam \\
        --num-threads ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
