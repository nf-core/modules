process GENRICH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genrich:0.6.1--h5bf99c6_1' :
        'biocontainers/genrich:0.6.1--h5bf99c6_1' }"

    input:
    tuple val(meta), path(treatment_bam), path(control_bam)
    path  blacklist_bed

    output:
    tuple val(meta), path("*.narrowPeak")                     , emit: peak
    path "versions.yml"                                       , emit: versions

    tuple val(meta), path("*.pvalues.bedGraph"), optional:true, emit: bedgraph_pvalues
    tuple val(meta), path("*.pileup.bedGraph") , optional:true, emit: bedgraph_pileup
    tuple val(meta), path("*.intervals.bed")   , optional:true, emit: bed_intervals
    tuple val(meta), path("*.duplicates.txt")  , optional:true, emit: duplicates

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ""
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def treatment  = treatment_bam   ? "-t ${treatment_bam.join(',')}"   : ""
    def control    = control_bam     ? "-c ${control_bam.join(',')}"     : ""
    def blacklist  = blacklist_bed   ? "-E $blacklist_bed"               : ""

    if (meta.single_end && (!args.contains("-y") && !args.contains("-w"))) {
        log.info '[Genrich] Single-end data can only be analyzed if unpaired alignments are kept (-y or -w <int>), defaulting to -y option.'
        args = "-y ${args}"
    }
    if (args.contains("-R") && !args.contains("-r")) {
        log.info '[Genrich] Duplicates can only be saved if they are filtered out, defaulting to -r option (Remove PCR duplicates).'
        args = "-r ${args}"
    }
    """
    Genrich \\
        $args \\
        $treatment \\
        $control \\
        $blacklist \\
        -o ${prefix}.narrowPeak

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genrich: \$(echo \$(Genrich --version 2>&1) | sed 's/^Genrich, version //; s/ .*\$//')
    END_VERSIONS
    """
}
