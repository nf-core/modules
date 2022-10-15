process MSISENSOR_SCAN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::msisensor=0.5" : null)
    def container_image = "msisensor:0.5--hb3646a4_2"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tab"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    msisensor \\
        scan \\
        -d $fasta \\
        -o ${prefix}.msisensor_scan.tab \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msisensor: \$(msisensor 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
