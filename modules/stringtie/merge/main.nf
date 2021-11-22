process STRINGTIE_MERGE {
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.1.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.7--h978d192_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.7--h978d192_0"
    }

    input:
    path stringtie_gtf
    path annotation_gtf

    output:
    path "stringtie.merged.gtf", emit: gtf
    path  "versions.yml"       , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    stringtie \\
        --merge $stringtie_gtf \\
        -G $annotation_gtf \\
        -o stringtie.merged.gtf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
