
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: discuss if how to add samtools
    conda (params.enable_conda ? "bioconda::pretextmap=0.1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextmap%3A0.1.9--h9f5acd7_1':
        'quay.io/biocontainers/pretextmap:0.1.9--h9f5acd7_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pretext"), emit: pretext
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
	-h \\
	$bam | PretextMap \\
        $args \\
        -o ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: PretextMap | grep "Version" | sed 's/PretextMap Version //g'
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
