process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ivar=1.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.3.1--h089eab3_0' :
        'quay.io/biocontainers/ivar:1.3.1--h089eab3_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta

    output:
    tuple val(meta), path("*.fa")      , emit: fasta
    tuple val(meta), path("*.qual.txt"), emit: qual
    tuple val(meta), path("*.mpileup") , optional:true, emit: mpileup
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def save_mpileup = params.save_mpileup ? "tee ${prefix}.mpileup |" : ""
    """
    samtools mpileup \\
        --reference $fasta \\
        $args2 \\
        $bam | \\
        $save_mpileup \\
        ivar consensus \\
            $args \\
            -p $prefix

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//')
    END_VERSIONS
    """
}
