process SNPDISTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snp-dists=0.8.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snp-dists:0.8.2--h5bf99c6_0"
    } else {
        container "quay.io/biocontainers/snp-dists:0.8.2--h5bf99c6_0"
    }

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    snp-dists \\
        $args \\
        $alignment > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
    """
}
