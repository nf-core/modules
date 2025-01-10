process ARRIBA_VISUALISATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::arriba=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.4.0--h6b7c446_1':
        'biocontainers/arriba:2.4.0--h6b7c446_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fusions)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(protein_domains)
    tuple val(meta4), path(cytobands)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cytobands = cytobands ? " --cytobands=$cytobands" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def protein_domains = protein_domains ? "--proteinDomains=$protein_domains" : ""
    """
    draw_fusions.R \\
        --fusions=$fusions \\
        --alignments=$bam \\
        --output=${prefix}.pdf \\
        --annotation=${gtf} \\
        $cytobands \\
        $protein_domains \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
