process DELLY_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::delly=0.8.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delly:0.8.7--he03298f_1' :
        'quay.io/biocontainers/delly:0.8.7--he03298f_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val(meta), path("*.csi"), emit: csi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    delly \\
        call \\
        $args \\
        -o ${prefix}.bcf \\
        -g  $fasta \\
        $bam \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}
