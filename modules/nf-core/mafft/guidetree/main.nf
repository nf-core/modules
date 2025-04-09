process MAFFT_GUIDETREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1':
        'biocontainers/mafft:7.525--h031d066_1' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mafft \\
        --thread ${task.cpus} \\
        --retree 0 \\
        ${args} \\
        --treeout ${fasta} > log.txt

    mv *.tree ${prefix}.dnd.tmp

    # remove all prefixes added by mafft which make the output incompatible with other tools
    awk '{gsub(/^[0-9]+_/, ""); print}' ${prefix}.dnd.tmp > ${prefix}.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
