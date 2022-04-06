process SAMTOOLS_AMPLICONCLIP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15--h1170115_1' :
        'quay.io/biocontainers/samtools:1.15--h1170115_1' }"

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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rejects = save_cliprejects ? "--rejects-file ${prefix}.cliprejects.bam" : ""
    def stats   = save_clipstats   ? "-f ${prefix}.clipstats.txt"               : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
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
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
