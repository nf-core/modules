process GFFCOMPARE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
        'biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"

    input:
    tuple val(meta), path(gtfs)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(reference_gtf)

    output:
    tuple val(meta), path("*.annotated.gtf"), optional: true, emit: annotated_gtf
    tuple val(meta), path("*.combined.gtf") , optional: true, emit: combined_gtf
    tuple val(meta), path("*.tmap")         , optional: true, emit: tmap
    tuple val(meta), path("*.refmap")       , optional: true, emit: refmap
    tuple val(meta), path("*.loci")         , emit: loci
    tuple val(meta), path("*.stats")        , emit: stats
    tuple val(meta), path("*.tracking")     , emit: tracking
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_fasta = fasta ? "-s ${fasta}" : ''
    def ref_gtf = reference_gtf ? "-r ${reference_gtf}" : ''
    """
    gffcompare \\
        $args \\
        $ref_fasta \\
        $ref_gtf \\
        -o $prefix \\
        $gtfs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(echo \$(gffcompare --version 2>&1) | sed 's/^gffcompare v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotated.gtf
    touch ${prefix}.combined.gtf
    touch ${prefix}.tmap
    touch ${prefix}.refmap
    touch ${prefix}.loci
    touch ${prefix}.stats
    touch ${prefix}.tracking

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(echo \$(gffcompare --version 2>&1) | sed 's/^gffcompare v//')
    END_VERSIONS
    """
}
