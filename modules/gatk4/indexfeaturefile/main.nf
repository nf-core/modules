process GATK4_INDEXFEATUREFILE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0' :
        'quay.io/biocontainers/gatk4:4.2.0.0--0' }"

    input:
    tuple val(meta), path(feature_file)

    output:
    tuple val(meta), path("*.{tbi,idx}"), emit: index
    path  "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    gatk \\
        IndexFeatureFile \\
        $args \\
        -I $feature_file

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
