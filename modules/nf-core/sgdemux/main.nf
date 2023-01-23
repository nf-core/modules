process SGDEMUX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::sgdemux=1.1.0"
    if (workflow.profile.tokenize(',').intersect(['Singularity', 'Docker']).size() >= 1) {
        exit 1, "SGDEMUX module does not support Singularity or Docker. Please use Conda / Mamba instead."
    }
    // TODO check sgdemux is not on singularity/docker (google says no) double check with Nils
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        // 'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        // 'quay.io/biocontainers/YOUR-TOOL-HERE' }"]

    input:
    tuple val(meta), path(run_manifest), path(run_dir)

    output:
    // TODO output/Samples/*/*_R*...
    tuple val(meta), path('output/*_R*.fastq.gz'), emit: sample_fastq
    tuple val(meta), path('output/metrics.tsv'), emit: metrics
    tuple val(meta), path('output/most_frequent_unmatched.tsv'), emit: most_frequent_unmatched
    tuple val(meta), path('output/per_project_metrics.tsv'), emit: per_project_metrics
    tuple val(meta), path('output/per_sample_metrics.tsv'), emit: per_sample_metrics
    tuple val(meta), path('output/sample_barcode_hop_metrics.tsv'), emit: sample_barcode_hop_metrics
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def runManifest = run_manifest ? "-r ${run_manifest}" : ""
    """
    mkdir -p output/
    sgdemux \\
        --sample-metadata $run_manifest \\
        --fastqs $run_dir \\
        --output-dir output/ \\
        --demux-threads $task.cpus \\
        --compressor-threads $task.cpus \\
        --writer-threads $task.cpus \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sgdemux: \$(echo \$(sgdemux --version 2>&1) | cut -d " " -f2)
    END_VERSIONS
    """
}