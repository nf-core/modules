process ILASTIK_MULTICUT {
    tag "$meta.id"
    label 'process_low'

    container "labsyspharm/mcmicro-ilastik:1.6.0"

    input:
    tuple val(meta), path(h5)
    path ilp
    path probs

    output:
    tuple val(meta), path("*.h5") , emit: h5
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
        --raw_data=$h5 \\
        --probabilities=$probs \\
        --export_source="Multicut Segmentation" \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(/ilastik-release/run_ilastik.sh --headless --version)
    END_VERSIONS
    """
}
