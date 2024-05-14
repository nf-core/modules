process BASES2FASTQ {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/bases2fastq:1.1.0"

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BASES2FASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def runManifest = run_manifest ? "-r ${run_manifest}" : ""
    """
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
