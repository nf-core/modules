process CELLRANGER_MKGTF {
    tag "$gtf"
    label 'process_low'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "nfcore/cellranger:6.1.2"

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
    cellranger \\
        mkgtf \\
        $gtf \\
        ${gtf.baseName}.filtered.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
