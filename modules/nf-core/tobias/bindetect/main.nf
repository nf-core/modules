process TOBIAS_BINDETECT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.17.5--py310h3479294_0':
        'quay.io/biocontainers/tobias:0.17.5--py310h3479294_0' }"

    input:
    tuple val(meta), path(signals), path(peaks), path(motifs), path(fasta)

    output:
    tuple val(meta), path("bindetect")                , emit: outdir
    tuple val(meta), path("bindetect/*_results.txt")  , emit: results, optional: true
    tuple val(meta), path("bindetect/*_results.xlsx") , emit: results_xlsx, optional: true
    tuple val(meta), path("bindetect/*_distances.txt"), emit: distances, optional: true
    tuple val(meta), path("bindetect/*_figures.pdf")  , emit: figures, optional: true
    tuple val("${task.process}"), val('tobias'), eval('TOBIAS --version'), topic: versions, emit: versions_tobias

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p .matplotlib
    export MPLCONFIGDIR="\${PWD}/.matplotlib"

    TOBIAS BINDetect \\
        --signals $signals \\
        --motifs $motifs \\
        --genome $fasta \\
        --peaks $peaks \\
        --outdir bindetect \\
        --prefix $prefix \\
        --cores ${task.cpus} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bindetect/TF1/beds
    touch bindetect/${prefix}_results.txt
    touch bindetect/${prefix}_results.xlsx
    touch bindetect/${prefix}_distances.txt
    touch bindetect/${prefix}_figures.pdf
    touch bindetect/TF1/TF1_overview.txt
    touch bindetect/TF1/beds/TF1_all.bed
    """
}
