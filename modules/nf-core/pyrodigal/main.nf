process PYRODIGAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2fe9a8ce513c91df34b43a6610df94c3a2eb3bd0:da1134ad604a59a6f439bdcc3f6df690eba47e9a-0':
        'biocontainers/mulled-v2-2fe9a8ce513c91df34b43a6610df94c3a2eb3bd0:da1134ad604a59a6f439bdcc3f6df690eba47e9a-0' }"

    input:
    tuple val(meta), path(fasta)
    val(output_format)

    output:
    tuple val(meta), path("*.${output_format}.gz")      , emit: annotations
    tuple val(meta), path("*.fna.gz")                   , emit: fna
    tuple val(meta), path("*.faa.gz")                   , emit: faa
    tuple val(meta), path("*.score.gz")                 , emit: score
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pigz -cdf ${fasta} > pigz_fasta.fna

    pyrodigal \\
        -j ${task.cpus} \\
        $args \\
        -i pigz_fasta.fna \\
        -f $output_format \\
        -o "${prefix}.${output_format}" \\
        -d ${prefix}.fna \\
        -a ${prefix}.faa \\
        -s ${prefix}.score

    pigz -nmf ${prefix}*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: \$(echo \$(pyrodigal --version 2>&1 | sed 's/pyrodigal v//'))
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${output_format}.gz
    touch ${prefix}.fna.gz
    touch ${prefix}.faa.gz
    touch ${prefix}.score.gz
    touch versions.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: \$(echo \$(pyrodigal --version 2>&1 | sed 's/pyrodigal v//'))
    END_VERSIONS
    """
}
