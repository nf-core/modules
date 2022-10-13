process FASTTREE {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fasttree=2.1.10" : null)
        def container_image = "/fasttree:2.1.10--h516909a_4"
                                           container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    path alignment

    output:
    path "*.tre",         emit: phylogeny
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fasttree \\
        $args \\
        -log fasttree_phylogeny.tre.log \\
        -nt $alignment \\
        > fasttree_phylogeny.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -1  | sed 's/^FastTree \\([0-9\\.]*\\) .*\$/\\1/')
    END_VERSIONS
    """
}
