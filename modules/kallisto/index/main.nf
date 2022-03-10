process KALLISTO_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1' :
        'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1' }"

    input:
    path fasta

    output:
    path "kallisto" , emit: idx
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto \\
        index \\
        $args \\
        -i kallisto \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """
}
