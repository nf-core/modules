process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'quay.io/biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    path gff

    output:
    path "*.gtf"        , emit: gtf
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${gff.baseName}${task.ext.suffix}" : "${gff.baseName}"
    """
    gffread \\
        $gff \\
        $args \\
        -o ${prefix}.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
