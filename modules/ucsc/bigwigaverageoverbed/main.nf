def VERSION = '377'

process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-bigwigaverageoverbed=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
    } else {
        container "quay.io/biocontainers/ucsc-bigwigaverageoverbed:377--h0b8a92a_2"
    }

    input:
    tuple val(meta), path(bed)
    path bigwig

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    # there is a bug that bigWigAverageOverBed can not handle ensembl seqlevels style.
    bigWigAverageOverBed \\
        $options.args \\
        $bigwig \\
        $bed \\
        ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
