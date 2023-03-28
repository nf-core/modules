process CELLRANGER_ARC_MKGTF {
    tag "$gtf"
    label 'process_low'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_MKGTF module does not support Conda.
        Please use docker or singularity containers."
    }
    container "heylf/cellranger-arc:2.0.2"

    input:
    path gtf

    output:
    path "*.filtered.gtf", emit: gtf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cellranger-arc \\
        mkgtf \\
        $gtf \\
        ${gtf.baseName}.filtered.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
