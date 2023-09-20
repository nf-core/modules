
process LEARNMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::learnmsa=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/learnmsa:1.3.0--pyhdfd78af_0':
        'biocontainers/learnmsa:1.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(fasta)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    learnMSA \\
        $args \\
        -i $fasta \\
        -o ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        learnMSA: \$(learnMSA -h | grep -oP 'version \\K[0-9.]+')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        learnMSA: \$(learnMSA -h | grep -oP 'version \\K[0-9.]+')
    END_VERSIONS
    """
}
