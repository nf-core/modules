process GENRICH {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::genrich=0.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genrich:0.6.1--h5bf99c6_1' :
        'quay.io/biocontainers/genrich:0.6.1--h5bf99c6_1' }"

    input:
    tuple val(meta), path(treatment_bam)
    path  control_bam
    path  blacklist_bed
    val   save_pvalues
    val   save_pileup
    val   save_bed
    val   save_duplicates

    output:
    tuple val(meta), path("*narrowPeak")                     , emit: peaks
    tuple val(meta), path("*pvalues.bedGraph"), optional:true, emit: bedgraph_pvalues
    tuple val(meta), path("*pileup.bedGraph") , optional:true, emit: bedgraph_pileup
    tuple val(meta), path("*intervals.bed")   , optional:true, emit: bed_intervals
    tuple val(meta), path("*duplicates.txt")  , optional:true, emit: duplicates
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control    = control_bam    ? "-c $control_bam"               : ''
    def blacklist  = blacklist_bed  ? "-E $blacklist_bed"             : ""
    def pvalues    = save_pvalues   ? "-f ${prefix}.pvalues.bedGraph" : ""
    def pileup     = save_pileup    ? "-k ${prefix}.pileup.bedGraph"  : ""
    def bed        = save_bed       ? "-b ${prefix}.intervals.bed"    : ""
    def duplicates = ""
    if (save_duplicates) {
        if (args.contains('-r')) {
            duplicates = "-R ${prefix}.duplicates.txt"
        } else {
            log.info '[Genrich] Duplicates can only be saved if they are filtered, defaulting to -r option (Remove PCR duplicates).'
            duplicates = "-r -R ${prefix}.duplicates.txt"
        }
    }
    """
    Genrich \\
        -t $treatment_bam \\
        $args \\
        $control \\
        $blacklist \\
        -o ${prefix}.narrowPeak \\
        $pvalues \\
        $pileup \\
        $bed \\
        $duplicates \\
        $control

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genrich: \$(echo \$(Genrich --version 2>&1) | sed 's/^Genrich, version //; s/ .*\$//')
    END_VERSIONS
    """
}
