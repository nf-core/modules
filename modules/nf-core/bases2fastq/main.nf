process BASES2FASTQ {
    tag "$meta.id"
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using bases2fastq. Please use docker or singularity containers."
    }
    container "elembio/bases2fastq:1.1.0"

    input:
    tuple val(meta), path(run_manifest), path(run_dir)

    output:
    tuple val(meta), path('output/Samples/*/*_R*.fastq.gz'), emit: sample_fastq
    tuple val(meta), path('output/Samples/*/*.json')       , emit: sample_json
    tuple val(meta), path('output/*.html')                 , emit: qc_report
    tuple val(meta), path('output/RunStats.json')          , emit: run_stats
    tuple val(meta), path('output/RunManifest.json')       , emit: generated_run_manifest
    tuple val(meta), path('output/Metrics.csv')            , emit: metrics
    tuple val(meta), path('output/UnassignedSequences.csv'), emit: unassigned
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def runManifest = run_manifest ? "-r ${run_manifest}" : ""
    """
    ls

    bases2fastq \\
        -p $task.cpus \\
        $runManifest \\
        $args \\
        $run_dir \\
        output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bases2fastq: \$(bases2fastq --version | sed -e "s/bases2fastq version //g")
    END_VERSIONS
    """
}
