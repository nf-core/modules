process YAHA_ALIGN {
    tag "$meta_reads.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::yaha=0.1.83 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3d7bae608c81b0a9427122d41c428d8c9696105:b765610227f847a2565270e300496a62199e0e31-0':
        'quay.io/biocontainers/mulled-v2-f3d7bae608c81b0a9427122d41c428d8c9696105:b765610227f847a2565270e300496a62199e0e31-0' }"

    input:
    tuple val(meta_reads), path(reads)
    tuple val(meta_compressed_genome), path(compressed_genome)
    tuple val(meta_index), path(index)

    output:
    tuple val(meta_reads), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_reads.id}"
    def threads = task.cpus
    """
    yaha \\
        -t $threads \\
        -x $index \\
        -q $reads \\
        $args -oss stdout | samtools view -Sb - > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yaha: \$(echo \$(yaha 2>&1) | sed 's/^.*YAHA version //; s/Usage.*\$//' )
    END_VERSIONS
    """
}
