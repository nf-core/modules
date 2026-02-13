process BASES2FASTQ {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/elembio/bases2fastq:2.3.0"

    input:
    tuple val(meta), path(run_manifest), path(run_dir, stageAs: 'input_dir')

    output:
    tuple val(meta), path("${prefix}/Samples/**/*_R*.fastq.gz"), emit: sample_fastq
    tuple val(meta), path("${prefix}/Samples/**/*_stats.json") , emit: sample_json
    tuple val(meta), path("${prefix}/*_QC.html")               , emit: qc_report
    tuple val(meta), path("${prefix}/multiqc_report.html")     , emit: multiqc_report, optional: true
    tuple val(meta), path("${prefix}/RunStats.json")           , emit: run_stats
    tuple val(meta), path("${prefix}/RunManifest.json")        , emit: generated_run_manifest
    tuple val(meta), path("${prefix}/Metrics.csv")             , emit: metrics
    tuple val(meta), path("${prefix}/UnassignedSequences.csv") , emit: unassigned
    tuple val("${task.process}"), val('bases2fastq'), eval('bases2fastq --version | sed "s/.*version //;s/,.*//"'), emit: versions_bases2fastq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BASES2FASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def runManifest = run_manifest ? "-r ${run_manifest}" : ""
    """
    bases2fastq \\
        -p $task.cpus \\
        $runManifest \\
        $args \\
        $run_dir \\
        ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/Samples/DefaultSample

    touch ${prefix}/Metrics.csv
    touch ${prefix}/RunManifest.json
    touch ${prefix}/UnassignedSequences.csv
    echo | gzip > ${prefix}/Samples/DefaultSample/DefaultSample_R1.fastq.gz
    echo | gzip > ${prefix}/Samples/DefaultSample/DefaultSample_R2.fastq.gz
    touch ${prefix}/Bases2Fastq-Sim_QC.html
    touch ${prefix}/multiqc_report.html
    touch ${prefix}/RunStats.json
    touch ${prefix}/Samples/DefaultSample/DefaultSample_stats.json
    """
}
