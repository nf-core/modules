// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf3ee0025f72df3ca5cabfae265b3623dcd9f9239f993ea258bf5e7846e9a4e4/data'
        : 'community.wave.seqera.io/library/ivar:1.4.4--60def7bdef0d4d2f'}"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fai
    path  gff
    val   save_mpileup

    output:
    tuple val(meta), path("*.tsv")    , emit: tsv
    tuple val(meta), path("*.mpileup"), optional:true, emit: mpileup
    tuple val("${task.process}"), val('ivar'), eval("ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p'"), emit: versions_ivar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def features = gff ? "-g $gff" : ""
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    """
    samtools \\
        mpileup \\
        $args2 \\
        --reference $fasta \\
        $bam \\
        $mpileup \\
        | ivar \\
            variants \\
            $args \\
            $features \\
            -r $fasta \\
            -p $prefix
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_mpileup = save_mpileup ? "touch ${prefix}.mpileup" : ''
    """
    touch ${prefix}.tsv
    $touch_mpileup
    """
}
