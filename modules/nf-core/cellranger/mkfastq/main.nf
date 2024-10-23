process CELLRANGER_MKFASTQ {
    tag {"$meta.lane" ? "$meta.id"+"."+"$meta.lane" : "$meta.id" }
    label 'process_medium'

    container "nf-core/cellrangermkfastq:8.0.0"

    input:
    tuple val(meta), path(csv), path(bcl)

    output:
    tuple val(meta), path("*_outs/outs/fastq_path/**/*.fastq.gz")          , emit: fastq
    tuple val(meta), path("*_outs/outs/fastq_path/Undetermined*.fastq.gz") , optional:true, emit: undetermined_fastq
    tuple val(meta), path("*_outs/outs/fastq_path/Reports")                , emit: reports
    tuple val(meta), path("*_outs/outs/fastq_path/Stats")                  , emit: stats
    tuple val(meta), path("*_outs/outs/interop_path/*.bin")                , emit: interop
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_MKFASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" //run_dir (bcl) and id must be different because a folder is created with the id value
    """
    cellranger \\
        mkfastq \\
        --id=${prefix}_outs \\
        --run=$bcl \\
        --csv=$csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # data for reports output channel
    mkdir -p "${prefix}_outs/outs/fastq_path/Reports"
    # data for stats output channel
    mkdir -p "${prefix}_outs/outs/fastq_path/Stats"
    # data for interops output channel
    mkdir -p "${prefix}_outs/outs/interop_path/"
    touch "${prefix}_outs/outs/interop_path/IndexMetricsOut.bin"
    # data for fastq channels
    mkdir -p "${prefix}_outs/outs/fastq_path/sample/files/"
    touch "${prefix}_outs/outs/fastq_path/Undetermined_fake_file.fastq"
    touch "${prefix}_outs/outs/fastq_path/sample/files/fake_file.fastq"
    gzip "${prefix}_outs/outs/fastq_path/Undetermined_fake_file.fastq"
    gzip "${prefix}_outs/outs/fastq_path/sample/files/fake_file.fastq"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
