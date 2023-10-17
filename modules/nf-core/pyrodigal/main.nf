process PYRODIGAL {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyrodigal=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2fe9a8ce513c91df34b43a6610df94c3a2eb3bd0:697b3838b186fac6a9ceec198b09d4032162a079-0':
        'biocontainers/mulled-v2-2fe9a8ce513c91df34b43a6610df94c3a2eb3bd0:697b3838b186fac6a9ceec198b09d4032162a079-0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gff.gz")      , emit: gff
    tuple val(meta), path("*.fna.gz")      , emit: fna
    tuple val(meta), path("*.faa.gz")      , emit: faa
    tuple val(meta), path("*.score.gz")    , emit: score
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pigz -cdf ${fasta} > pigz_fasta.fna

    pyrodigal \\
        $args \\
        -i pigz_fasta.fna \\
        -o ${prefix}.gff \\
        -d ${prefix}.fna \\
        -a ${prefix}.faa \\
        -s ${prefix}.score

    pigz -nmf ${prefix}*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: \$(echo \$(pyrodigal --version 2>&1 | sed 's/pyrodigal v//'))
    END_VERSIONS
    """
}
