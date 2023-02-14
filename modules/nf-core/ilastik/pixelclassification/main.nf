process ILASTIK_PIXELCLASSIFICATION {
    tag "$meta.id"
    label 'process_single'

    container "labsyspharm/mcmicro-ilastik:1.5.1"

    input:
    tuple val(meta), path(h5)
    path ilp

    output:
    tuple val(meta), path("*.h5"), emit: h5
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /ilastik-release/run_ilastik.sh \\
        --headless \\
        --readonly 1 \\
        --project=$ilp \\
        $args \\
        $h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(/ilastik-release/run_ilastik.sh --headless --version)
    END_VERSIONS
    """
}
