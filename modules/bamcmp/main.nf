def VERSION = '2.2'
process BAMCMP {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bamcmp=2.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamcmp:2.2--h05f6578_0' :
        'quay.io/biocontainers/bamcmp:2.2--h05f6578_0' }"

    input:
    tuple val(meta), path(sample), path(contaminant)

    output:
    tuple val(meta), path("*primary.bam")      , emit: bam
    tuple val(meta), path("*contamination.bam"), emit: contamination_bam
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamcmp \\
        -s "as" \\
        -1 $sample \\
        -2 $contaminant \\
        -A ${prefix}_primary.bam \\
        -B ${prefix}_contamination.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamcmp: $VERSION
    END_VERSIONS
    """

}
