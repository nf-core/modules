process HYPO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hypo=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hypo:1.0.3--h9a82719_1':
        'biocontainers/hypo:1.0.3--h9a82719_1' }"

    input:
    tuple val(meta), path(sr_bam)
    tuple val(meta2), path(reads)
    path draft
    val genome_size
    val reads_coverage

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "${reads.join(' ')}" | tr " " "\\n" > sr.fofn
    hypo \\
        -r @sr.fofn \\
        -d $draft \\
        -b $sr_bam \\
        -c $reads_coverage \\
        -s $genome_size \\
        -t $task.cpus \\
        -o ${prefix}.fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hypo: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: See above
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hypo: $VERSION
    END_VERSIONS
    """
}
