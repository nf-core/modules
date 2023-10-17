process MAFFT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::mafft=7.520"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.520--hec16e2b_1':
        'biocontainers/mafft:7.520--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)
    path  addsequences

    output:
    tuple val(meta), path("*.fas"), emit: fas
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def add = addsequences ? "--add $addsequences" : ''
    """
    mafft \\
        --thread ${task.cpus} \\
        ${args} \\
        ${add} \\
        ${fasta} \\
        > ${prefix}.fas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
