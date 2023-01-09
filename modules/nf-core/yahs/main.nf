process YAHS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::yahs=1.2a.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yahs:1.2a.2--h7132678_0':
        'quay.io/biocontainers/yahs:1.2a.2--h7132678_0' }"

    input:
    tuple val(meta), path(hic_map)
    path fasta
    path fai

    output:
    tuple val(meta), path("*scaffolds_final.fa") , emit: scaffolds_fasta
    tuple val(meta), path("*scaffolds_final.agp"), emit: scaffolds_agp
    tuple val(meta), path("*bin")                , emit: binary
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    yahs $args \\
        -o $prefix \\
        $fasta \\
        $hic_map

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version 2>&1)
    END_VERSIONS
    """
}
