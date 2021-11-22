process EMMTYPER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::emmtyper=0.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/emmtyper:0.2.0--py_0"
    } else {
        container "quay.io/biocontainers/emmtyper:0.2.0--py_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    emmtyper \\
        $options.args \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' )
    END_VERSIONS
    """
}
