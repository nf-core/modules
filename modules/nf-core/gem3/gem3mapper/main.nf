process GEM3_GEM3MAPPER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gem3-mapper=3.6.1 bioconda::samtools=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-240a9c1936dd6a68f46aa198b2629b6734a18428:8fc2f8fd66e104a88e3eec4d8fec1665ee7f7278-0':
        'biocontainers/mulled-v2-240a9c1936dd6a68f46aa198b2629b6734a18428:8fc2f8fd66e104a88e3eec4d8fec1665ee7f7278-0' }"

    input:
    tuple val(meta), path(index)
    tuple val(meta2), path(fastq)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gem-mapper -F 'SAM' -I $index -i $fastq -t $task.cpus | samtools view -bS -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem-mapper: \$(echo \$(gem-mapper --version 2>&1) | sed 's/v//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem-mapper: \$(echo \$(gem-mapper --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}