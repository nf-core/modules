process GLNEXUS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::glnexus=1.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glnexus:1.4.1--h40d77a6_0' :
        'quay.io/biocontainers/glnexus:1.4.1--h40d77a6_0' }"

    input:
    tuple val(meta), path(gvcfs)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

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
    "${task.process}":
        glnexus: \$( echo \$(glnexus_cli 2>&1) | head -n 1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$( echo \$(glnexus_cli 2>&1) | head -n 1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """
}
