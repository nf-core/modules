process COOLER_DIGEST {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    def container_image = "cooler:0.8.11--pyh3252c3a_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    path fasta
    path chromsizes
    val  enzyme

    output:
    path "*.bed"                  , emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cooler digest \\
        $args \\
        -o "${fasta.baseName}_${enzyme.replaceAll(/[^0-9a-zA-Z]+/, '_')}.bed" \\
        $chromsizes \\
        $fasta \\
        $enzyme

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
