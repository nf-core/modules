
process MMSEQS_EASYCLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    tuple val(meta), path(sequence)

    output:
    tuple val(meta), path("*rep_seq.fasta") , emit: representatives
    tuple val(meta), path("*all_seqs.fasta"), emit: fasta
    tuple val(meta), path("*.tsv")          , emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mmseqs \\
        easy-cluster \\
        ${sequence} \\
        ${prefix} \\
        tmp1 \\
        $args \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.tsv
    touch ${prefix}_rep_seq.fasta
    touch ${prefix}_all_seqs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
