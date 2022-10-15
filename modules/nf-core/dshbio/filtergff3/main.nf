process DSHBIO_FILTERGFF3 {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::dsh-bio=2.1" : null)
    def container_image = "dsh-bio:2.1--hdfd78af_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.gff3.gz"), emit: gff3
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dsh-bio \\
        filter-gff3 \\
        $args \\
        -i $gff3 \\
        -o ${prefix}.gff3.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dshbio: \$(dsh-bio --version 2>&1 | grep -o 'dsh-bio-tools .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
