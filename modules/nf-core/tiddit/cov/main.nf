process TIDDIT_COV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::tiddit=3.0.0" : null)
    def container_image = "/tiddit:3.0.0--py39h59fae87_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(input)
    path  fasta

    output:
    tuple val(meta), path("*.bed"), optional: true, emit: cov
    tuple val(meta), path("*.wig"), optional: true, emit: wig
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref $fasta" : ""
    """
    tiddit \\
        --cov \\
        -o $prefix \\
        $args \\
        --bam $input \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.wig
    touch ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiddit: \$(echo \$(tiddit 2>&1) | sed 's/^.*tiddit-//; s/ .*\$//')
    END_VERSIONS
    """
}
