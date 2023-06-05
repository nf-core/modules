process GANGSTR {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gangstr=2.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gangstr:2.5.0--h48cf4b7_4':
        'biocontainers/gangstr:2.5.0--h48cf4b7_4' }"

    input:
    tuple val(meta), path(alignment_files), path(alignment_indices), path(ref_regions)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf")              , emit: vcf
    tuple val(meta), path("*.samplestats.tab")  , emit: samplestats
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input = alignment_files.join(",")

    """
    GangSTR \\
        --bam ${input} \\
        --ref ${fasta} \\
        --regions ${ref_regions} \\
        --out ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gangstr: \$(echo \$(GangSTR --version 2>&1))
    END_VERSIONS
    """
}
