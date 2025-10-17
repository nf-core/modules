process ARTIC_GUPPYPLEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5a747cc579edfc0cb2176b749afc02550ab5de678ae6a40d2cfadeba6c0de25d/data' :
        'community.wave.seqera.io/library/artic:1.6.2--d4956cdc155b8612' }"

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    artic \\
        guppyplex \\
        ${args} \\
        --directory ${fastq_dir} \\
        --output ${prefix}.fastq

    pigz -p ${task.cpus} *.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
