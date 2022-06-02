process SKESA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::skesa=2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skesa:2.4.0--he1c1bb9_0':
        'quay.io/biocontainers/skesa:2.4.0--he1c1bb9_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contigs.fasta.gz')      , emit: contigs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory
    def illumina_reads = meta.single_end ? "${reads}" : "${reads[0]},${reads[1]}"
    if (meta.single_end == false) {
        args += " --use_paired_ends"
    }
    """
    skesa \\
    $args \\
    --cores $task.cpus \\
    --memory $maxmem \\
    --reads $illumina_reads \\
    | gzip > ${prefix}.contigs.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skesa: \$(skesa --version 2>&1 | tail -1 | awk '{ print \$2 }' | sed 's/v\\.//')
    END_VERSIONS
    """
}
