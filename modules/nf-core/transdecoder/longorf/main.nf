process TRANSDECODER_LONGORF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::transdecoder=5.5.0" : null)
    def container_image = "/transdecoder:5.5.0--pl5262hdfd78af_4"
    container { (params.container_registry ?: 'quay.io/comp-bio-aging' + container_image) }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*.pep") , emit: pep
    tuple val(meta), path("${meta.id}/*.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}/*.cds") , emit: cds
    tuple val(meta), path("${meta.id}/*.dat") , emit: dat
    path("${meta.id}/")                       , emit: folder
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    TransDecoder.LongOrfs \\
        $args \\
        -O $prefix \\
        -t \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.LongOrfs --version) | sed -e "s/TransDecoder.LongOrfs //g")
    END_VERSIONS
    """
}
