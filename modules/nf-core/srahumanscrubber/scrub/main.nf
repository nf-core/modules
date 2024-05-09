process SRAHUMANSCRUBBER_SCRUB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.0.0--hdfd78af_0':
        'quay.io/biocontainers/sra-human-scrubber:2.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), emit: scrubbed_reads
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if (meta.single_end) {
        """
        zcat ${reads} | scrub.sh -d $db | gzip > ${prefix}.scrubbed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: $VERSION
        END_VERSIONS
        """
    } else {
        """
        zcat ${reads[0]} | scrub.sh -d $db | gzip > ${prefix}_R1.scrubbed.fastq.gz
        zcat ${reads[1]} | scrub.sh -d $db | gzip > ${prefix}_R2.scrubbed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: $VERSION
            sra-human-scrubber-db: \$DBVERSION
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.0.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "" | zcat > "${prefix}.scrubbed.fastq.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: $VERSION
    END_VERSIONS
    """
}
