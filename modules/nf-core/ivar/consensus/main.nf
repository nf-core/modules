process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.4.3--h43eeafb_0' :
        'biocontainers/ivar:1.4.3--h43eeafb_0' }"

    input:
    tuple val(meta), path(bam)
    path fasta
    val save_mpileup

    output:
    tuple val(meta), path("*.fa")      , emit: fasta
    tuple val(meta), path("*.qual.txt"), emit: qual
    tuple val(meta), path("*.mpileup") , optional:true, emit: mpileup
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    """
    samtools \\
        mpileup \\
        --reference $fasta \\
        $args2 \\
        $bam \\
        $mpileup \\
        | ivar \\
            consensus \\
            $args \\
            -p $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_mpileup = save_mpileup ? "touch ${prefix}.mpileup" : ''
    """
    touch ${prefix}.fa
    touch ${prefix}.qual.txt
    $touch_mpileup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p')
    END_VERSIONS
    """
}
