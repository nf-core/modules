process YAHS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yahs:1.2--he4a0461_1':
        'biocontainers/yahs:1.2--he4a0461_1' }"

    input:
    tuple val(meta), path(hic_map)
    path fasta
    path fai

    output:
    tuple val(meta), path("*scaffolds_final.fa") , emit: scaffolds_fasta,  optional: true
    tuple val(meta), path("*scaffolds_final.agp"), emit: scaffolds_agp,    optional: true
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

    stub:
    """
    touch ${prefix}_scaffold_final.fa
    touch ${prefix}_scaffolds_final.agp
    touch ${prefix}.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version 2>&1)
    END_VERSIONS
    """
}
