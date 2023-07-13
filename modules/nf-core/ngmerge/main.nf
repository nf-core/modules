process NGMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ngmerge=0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngmerge:0.3--ha92aebf_1':
        'biocontainers/ngmerge:0.3--ha92aebf_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fq.gz"), emit: merged_reads
    tuple val(meta), path("*_1.fastq.gz")  , emit: unstitched_read1
    tuple val(meta), path("*_2.fastq.gz")  , emit: unstitched_read2
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    NGmerge \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -o ${prefix}.merged.fq.gz \\
        -f ${prefix}_unstitched \\
        -z \\
        -n $task.cpus \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGmerge: \$(echo \$(NGmerge --version 2>&1) | sed 's/^.*NGmerge, version //; s/ Copyright.*// ; s/: //g' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.merged.fq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGmerge: \$(echo \$(NGmerge --version 2>&1) | sed 's/^.*NGmerge, version //; s/ Copyright.*// ; s/: //g' ))
    END_VERSIONS
    """
}
