process GATK4_LEARNREADORIENTATIONMODEL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("*.tar.gz"), emit: artifactprior
    path "versions.yml"              , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def inputs_list = []
    f1r2.each() { a -> inputs_list.add(" -I " + a) }
    """
    gatk \\
        LearnReadOrientationModel \\
        ${inputs_list.join(' ')} \\
        -O ${prefix}.tar.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
