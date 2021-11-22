process GATK4_REVERTSAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gatk RevertSam \\
        I=$bam \\
        O=${prefix}.reverted.bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
