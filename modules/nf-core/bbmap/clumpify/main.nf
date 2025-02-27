process BBMAP_CLUMPIFY {
    tag "$meta.id"
    label 'process_single'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.17--h44cf5ac_0':
        'biocontainers/bbmap:39.17--he5f24ec_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=$reads" : "in1=${reads[0]} in2=${reads[1]}"
    def clumped  = meta.single_end ? "out=${prefix}.clumped.fastq.gz" : "out1=${prefix}_1.clumped.fastq.gz out2=${prefix}_2.clumped.fastq.gz"
    """
    clumpify.sh \\
        $raw \\
        $clumped \\
        $args \\
        &> ${prefix}.clumpify.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
