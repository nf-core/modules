process HIFIADAPTERFILT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiadapterfilt:3.0.0--hdfd78af_0':
        'biocontainers/hifiadapterfilt:3.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.filt.fastq.gz")       , emit: filt
    tuple val(meta), path("*.contaminant.blastout"), emit: blast_search
    tuple val(meta), path("*.stats")               , emit: stats
    tuple val(meta), path("*.blocklist")           , emit: headers
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The tool AUTOMATICALLY detects fastq files from the input folder, hence an explicit call of "fastq" is not needed.

    """
    hifiadapterfilt.sh \\
        ${args}

    mv *.filt.fastq.gz ${prefix}.filt.fastq.gz
    mv *.contaminant.blastout ${prefix}.contaminant.blastout
    mv *.stats ${prefix}.stats
    mv *.blocklist ${prefix}.blocklist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiadapterfilt: \$(hifiadapterfilt.sh -v)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.filt.fastq.gz
    touch ${prefix}.contaminant.blastout
    touch ${prefix}.stats
    touch ${prefix}.blocklist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiadapterfilt: \$(hifiadapterfilt.sh -v)
    END_VERSIONS
    """
}
