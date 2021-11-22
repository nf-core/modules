def VERSION = '377'

process UCSC_LIFTOVER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ucsc-liftover=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-liftover:377--h0b8a92a_3"
    } else {
        container "quay.io/biocontainers/ucsc-liftover:377--h0b8a92a_3"
    }

    input:
    tuple val(meta), path(bed)
    path(chain)

    output:
    tuple val(meta), path("*.lifted.bed")  , emit: lifted
    tuple val(meta), path("*.unlifted.bed"), emit: unlifted
    path "versions.yml"                    , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    liftOver \\
        $options.args \
        $bed \\
        $chain \\
        ${prefix}.lifted.bed \\
        ${prefix}.unlifted.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo "$VERSION")
    END_VERSIONS
    """
}
