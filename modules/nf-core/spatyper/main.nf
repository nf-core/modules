process SPATYPER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::spatyper=0.3.3" : null)
    def container_image = "/spatyper:0.3.3--pyhdfd78af_3"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fasta)
    path repeats
    path repeat_order

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_args = repeats && repeat_order ? "-r ${repeats} -o ${repeat_order}" : ""
    """
    spaTyper \\
        $args \\
        $input_args \\
        --fasta $fasta \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatyper: \$( echo \$(spaTyper --version 2>&1) | sed 's/^.*spaTyper //' )
    END_VERSIONS
    """
}
