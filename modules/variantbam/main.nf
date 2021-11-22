def VERSION = '1.4.4a'

process VARIANTBAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::variantbam=1.4.4a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/variantbam:1.4.4a--h7d7f7ad_5"
    } else {
        container "quay.io/biocontainers/variantbam:1.4.4a--h7d7f7ad_5"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam")         , emit: bam
    path "versions.yml"                    , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    variant \\
        $bam \\
        -o ${prefix}.bam \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
