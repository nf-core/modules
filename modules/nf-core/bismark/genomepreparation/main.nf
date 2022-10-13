process BISMARK_GENOMEPREPARATION {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
        def container_image = "/bismark:0.23.0--0"
                                                            container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path fasta, stageAs: "BismarkIndex/*"

    output:
    path "BismarkIndex" , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark_genome_preparation \\
        $args \\
        BismarkIndex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
