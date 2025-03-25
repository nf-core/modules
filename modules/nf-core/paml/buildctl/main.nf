process PAML_BUILDCTL {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(phy)
    path tree

    output:
    tuple val(meta), path("settings.ctl")   , emit: ctl
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "seqfile = ${phy}" > settings.ctl
    echo "treefile = ${tree}" >> settings.ctl
    echo "outfile = result.txt" >> settings.ctl
    echo "$args" >> settings.ctl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch settings.ctl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version)
    END_VERSIONS
    """
}
