process SAMTOOLS_AMPLICONCLIP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"

    input:
    tuple val(meta), path(bam)
    path bed
    val save_cliprejects
    val save_clipstats

    output:
    tuple val(meta), path("*.bam")            , emit: bam
    tuple val(meta), path("*.clipstats.txt")  , optional:true, emit: stats
    tuple val(meta), path("*.cliprejects.bam"), optional:true, emit: rejects_bam
    path "versions.yml"                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def rejects  = save_cliprejects ? "--rejects-file ${prefix}.cliprejects.bam" : ""
    def stats    = save_clipstats   ? "-f ${prefix}.clipstats.txt"               : ""
    """
    samtools \\
        ampliconclip \\
        $args \\
        $rejects \\
        $stats \\
        -b $bed \\
        -o ${prefix}.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
