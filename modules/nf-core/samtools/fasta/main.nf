process SAMTOOLS_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    def container_image = "samtools:1.15.1--h1170115_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def endedness = meta.single_end ? "-0 ${prefix}.fasta.gz" : "-1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz"
    """
    samtools \\
        fasta \\
        $args \\
        --threads ${task.cpus-1} \\
        $endedness \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
