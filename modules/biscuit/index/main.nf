process BISCUIT_INDEX {
    tag "$fasta"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.1.20211018" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.0.1.20211018--h81a5ba2_1':
        'quay.io/biocontainers/biscuit:1.0.1.20211018--h81a5ba2_1' }"

    input:
    path fasta, stageAs: "BiscuitIndex/*"

    output:
    path "BiscuitIndex", emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    biscuit \\
        index \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$(echo \$(biscuit version 2>&1) | sed 's/^.*BISCUIT Version: //; s/Using.*\$//')
    END_VERSIONS
    """
}
