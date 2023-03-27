process CELLRANGER_ARC_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_MKFASTQ module does not support Conda. 
        Please use docker or singularity containers."
    }
    container "nfcore/cellranger-arcmkfastq:2.0.2"

    input:
    path bcl
    path csv

    output:
    path "versions.yml", emit: versions
    path "${bcl.getSimpleName()}/outs/fastq_path/*.fastq.gz"  , emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cellranger-arc mkfastq --id=${bcl.getSimpleName()} \
        --run=$bcl \
        --csv=$csv \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${bcl.getSimpleName()}/outs/fastq_path/"
    touch ${bcl.getSimpleName()}/outs/fastq_path/fake_file.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
