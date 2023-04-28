process CELLRANGER_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    container "nfcore/cellrangermkfastq:7.1.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_MKFASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path bcl
    path csv

    output:
    path "**/outs/fastq_path/*.fastq.gz", emit: fastq
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bcl.getSimpleName()}"
    """
    cellranger \\
        mkfastq \\
        --id=${prefix} \\
        --run=$bcl \\
        --csv=$csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${bcl.getSimpleName()}"
    """
    mkdir -p "${prefix}/outs/fastq_path/"
    touch ${prefix}/outs/fastq_path/fake_file.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
