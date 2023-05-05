process PYRODIGAL {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pyrodigal=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyrodigal:2.1.0--py310h1425a21_0':
        'biocontainers/pyrodigal:2.1.0--py310h1425a21_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gff")      , emit: gff
    tuple val(meta), path("*.fna")      , emit: fna
    tuple val(meta), path("*.faa")      , emit: faa
    tuple val(meta), path("*.score")    , emit: score
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pyrodigal \\
        $args \\
        -i ${fasta} \\
        -o ${prefix}.gff \\
        -d ${prefix}.fna \\
        -a ${prefix}.faa \\
        -s ${prefix}.score

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: \$(echo \$(pyrodigal --version 2>&1 | sed 's/pyrodigal v//'))
    END_VERSIONS
    """
}
