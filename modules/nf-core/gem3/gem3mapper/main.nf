process GEM3_GEM3MAPPER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-240a9c1936dd6a68f46aa198b2629b6734a18428:543223d3cc2f69d86af72f7e9a3200812ae25327-0':
        'biocontainers/mulled-v2-240a9c1936dd6a68f46aa198b2629b6734a18428:543223d3cc2f69d86af72f7e9a3200812ae25327-0' }"

    input:
    tuple val(meta), path(gem)
    tuple val(meta2), path(fastq)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    gem-mapper -F 'SAM' -I $gem -i $fastq -t $task.cpus $args | samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem-mapper: \$(echo \$(gem-mapper --version 2>&1) | sed 's/v//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem-mapper: \$(echo \$(gem-mapper --version 2>&1) | sed 's/v//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
