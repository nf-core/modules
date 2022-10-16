process GENMAP_MAPPABILITY {
    tag "$index"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::genmap=1.3.0" : null)
    def container_image = "genmap:1.3.0--h1b792b2_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    input:
    path index

    output:
    path "*.wig"        , optional:true, emit: wig
    path "*.bedgraph"   , optional:true, emit: bedgraph
    path "*.txt"        , optional:true, emit: txt
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    genmap \\
        map \\
        $args \\
        -I $index \\
        -O mappability

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version 2>&1 | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """
}
