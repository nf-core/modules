process PRODIGAL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::prodigal=2.6.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prodigal:2.6.3--h516909a_2' :
        'quay.io/biocontainers/prodigal:2.6.3--h516909a_2' }"

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.${output_format}"), emit: gene_annotations
    tuple val(meta), path("${prefix}.fna"), emit: nucleotide_fasta
    tuple val(meta), path("${prefix}.faa"), emit: amino_acid_fasta
    tuple val(meta), path("${prefix}_all.txt"), emit: all_gene_annotations
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    prodigal -i "${genome}" \\
        $args \\
        -f $output_format \\
        -d "${prefix}.fna" \\
        -o "${prefix}.${output_format}" \\
        -a "${prefix}.faa" \\
        -s "${prefix}_all.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS
    """
}
