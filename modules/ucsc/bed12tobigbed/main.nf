def VERSION = '377'

process UCSC_BED12TOBIGBED {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-bedtobigbed=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:377--h446ed27_1"
    } else {
        container "quay.io/biocontainers/ucsc-bedtobigbed:377--h446ed27_1"
    }

    input:
    tuple val(meta), path(bed)
    path  sizes

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedToBigBed \\
        $bed \\
        $sizes \\
        ${prefix}.bigBed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
