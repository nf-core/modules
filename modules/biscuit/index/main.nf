process BISCUIT_INDEX {
    tag "$fasta"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::biscuit=1.0.2.20220113" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biscuit:1.0.2.20220113--h81a5ba2_0':
        'quay.io/biocontainers/biscuit:1.0.2.20220113--h81a5ba2_0' }"

    input:
    path fasta, stageAs: "BiscuitIndex/*"

    output:
    path "BiscuitIndex/*.fa*", emit: index, includeInputs: true
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    biscuit \\
        index \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
    END_VERSIONS
    """
}
