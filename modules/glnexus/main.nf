process GLNEXUS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::glnexus=1.4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/glnexus:1.4.1--h40d77a6_0"
    } else {
        container "quay.io/biocontainers/glnexus:1.4.1--h40d77a6_0"
    }

    input:
    tuple val(meta), path(gvcfs)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // Make list of GVCFs to merge
    def input = gvcfs.collect { it.toString() }
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Glnexus] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    glnexus_cli \\
        --threads $task.cpus \\
        --mem-gbytes $avail_mem \\
        $args \\
        ${input.join(' ')} \\
        > ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(glnexus_cli 2>&1) | head -n 1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """
}
