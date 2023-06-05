process CANU {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::canu=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'biocontainers/canu:2.2--ha47f30e_0' }"

    input:
    tuple val(meta), path(reads)
    val mode
    val genomesize

    output:
    tuple val(meta), path("*.report")                   , emit: report
    tuple val(meta), path("*.contigs.fasta.gz")         , emit: assembly                , optional: true
    tuple val(meta), path("*.unassembled.fasta.gz")     , emit: contigs
    tuple val(meta), path("*.correctedReads.fasta.gz")	, emit: corrected_reads         , optional: true
    tuple val(meta), path("*.trimmedReads.fasta.gz")	, emit: corrected_trimmed_reads , optional: true
    tuple val(meta), path("*.contigs.layout")           , emit: metadata                , optional: true
    tuple val(meta), path("*.contigs.layout.readToTig") , emit: contig_position         , optional: true
    tuple val(meta), path("*.contigs.layout.tigInfo")   , emit: contig_info             , optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["-pacbio", "-nanopore", "-pacbio-hifi"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Canu. Options: ${valid_mode.join(', ')}" }
    """
    canu \\
        -p ${prefix} \\
        $mode \\
        genomeSize=${genomesize} \\
        $args \\
        maxThreads=$task.cpus \\
        $reads

    gzip *.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """
}
