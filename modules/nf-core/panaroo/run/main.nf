process PANAROO_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::panaroo=1.2.9" : null)
    def container_image = "/panaroo:1.2.9--pyhdfd78af_0"
                                              container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                                      , emit: results
    tuple val(meta), path("results/core_gene_alignment.aln"), optional: true, emit: aln
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    panaroo \\
        $args \\
        -t $task.cpus \\
        -o results \\
        -i $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' ))
    END_VERSIONS
    """
}
