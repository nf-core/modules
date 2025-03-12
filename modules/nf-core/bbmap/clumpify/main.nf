process BBMAP_CLUMPIFY {
    tag "$meta.id"
    label 'process_single'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'wave.seqera.io/wt/5ce851926d93/wave/build:bbmap-39.18_pigz-2.8--4a8254d490b5baa0' }"

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
