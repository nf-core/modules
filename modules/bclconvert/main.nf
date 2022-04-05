process BCLCONVERT {
    tag 'demultiplexing'
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using bcl-convert. Please use docker or singularity containers."
    }
    container "nfcore/bclconvert:3.9.3"

    input:
    path samplesheet
    path run_dir

    output:
    path "*.fastq.gz"               ,emit: fastq
    path "Reports/*.{csv,xml,bin}"  ,emit: reports
    path "Logs/*.{log,txt}"         ,emit: logs
    path "versions.yml"             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcl-convert \
        $args \\
        --output-directory \$PWD \
        --bcl-input-directory ${run_dir} \\
        --sample-sheet ${samplesheet} \\
        --bcl-num-parallel-tiles $task.cpus \\
        --bcl-num-conversion-threads $task.cpus \\
        --bcl-num-compression-threads $task.cpus \\
        --bcl-num-decompression-threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bclconvert: \$(echo \$(bcl-convert -V 2>&1) | sed ''s/^.*Version //;s/Copyright.*//'' ))
    END_VERSIONS
    """
}
