process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::methyldackel=0.6.0' : null)
    def container_image = "/methyldackel:0.6.0--h22771d5_0"
                                                     container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.mbias.txt"), emit: txt
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MethylDackel mbias \\
        $args \\
        $fasta \\
        $bam \\
        $prefix \\
        --txt \\
        > ${prefix}.mbias.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methyldackel: \$(MethylDackel --version 2>&1 | cut -f1 -d" ")
    END_VERSIONS
    """
}
