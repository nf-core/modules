process CELLRANGER_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "nfcore/cellranger:6.1.2"

    input:
    path bcl
    path csv

    output:
    path "versions.yml", emit: versions
    path "${bcl.getSimpleName()}/outs/fastq_path/*.fastq.gz"  , emit: fastq

    script:
    def args = task.ext.args ?: ''
    """
    cellranger mkfastq --id=${bcl.getSimpleName()} \
        --run=$bcl \
        --csv=$csv
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
