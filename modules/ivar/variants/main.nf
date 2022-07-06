process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0' :
        'quay.io/biocontainers/ivar:1.3.1--h089eab3_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fai
    path  gff
    val   save_mpileup

    output:
    tuple val(meta), path("*.tsv")    , emit: tsv
    tuple val(meta), path("*.mpileup"), optional:true, emit: mpileup
    path "versions.yml"               , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//')
    END_VERSIONS
    """
}
