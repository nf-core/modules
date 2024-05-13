process CELLRANGER_MKFASTQ {
    tag "mkfastq"
    label 'process_medium'

    container "nf-core/cellrangermkfastq:8.0.0"

    input:
    path bcl
    path csv

    output:
    path "**/outs/fastq_path/*.fastq.gz", emit: fastq
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKFASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKFASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${bcl.getSimpleName()}"
    """
    mkdir -p "${prefix}/outs/fastq_path/"
    # data with something to avoid breaking nf-test java I/O stream
    cat <<-FAKE_FQ > ${prefix}/outs/fastq_path/fake_file.fastq
    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    +
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
    FAKE_FQ
    gzip -n ${prefix}/outs/fastq_path/fake_file.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
