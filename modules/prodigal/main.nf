process PRODIGAL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "prodigal=2.6.3 pigz=2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0' :
        'quay.io/biocontainers/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0' }"

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.${output_format}"),    emit: gene_annotations
    tuple val(meta), path("${prefix}.fna"),                 emit: nucleotide_fasta
    tuple val(meta), path("${prefix}.faa"),                 emit: amino_acid_fasta
    tuple val(meta), path("${prefix}_all.txt"),             emit: all_gene_annotations
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    pigz -cdf ${genome} | prodigal \\
        $args \\
        -f $output_format \\
        -d "${prefix}.fna" \\
        -o "${prefix}.${output_format}" \\
        -a "${prefix}.faa" \\
        -s "${prefix}_all.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
        pigz: \$(pigz -V 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """
}
