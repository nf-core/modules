process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kleborate=2.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kleborate:2.1.0--pyhdfd78af_1"
    } else {
        container "quay.io/biocontainers/kleborate:2.1.0--pyhdfd78af_1"
    }

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    kleborate \\
        $args \\
        --outfile ${prefix}.results.txt \\
        --assemblies *.fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
