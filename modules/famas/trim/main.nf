def VERSION = '0.0.7' // Version information provided by tool on CLI BUT I AM LAZY to parse it...so deal with it

process FAMAS_TRIM {
    tag "$meta.id"
    label 'process_low'

	// Nope...Docker only
    container 'solyris/famas:0.0.7'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    famas \\
        -i ${reads[0]} \\
        -j ${reads[1]} \\
        -o ${meta.id}.R1.trimmed.fastq.gz \\
        -p ${meta.id}.R2.trimmed.fastq.gz \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famas: $VERSION
    END_VERSIONS
    """
}