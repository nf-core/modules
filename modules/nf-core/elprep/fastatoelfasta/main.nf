process ELPREP_FASTATOELFASTA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/elprep:5.1.3--he881be0_1':
        'biocontainers/elprep:5.1.3--he881be0_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.elfasta")          , emit: elfasta
    tuple val(meta), path("logs/elprep/elprep*"), emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    elprep fasta-to-elfasta \\
        $fasta \\
        ${prefix}.elfasta \\
        --log-path ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(elprep 2>&1 | head -n2 | tail -n1 |sed 's/^.*version //;s/ compiled.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def timestamp = "${java.time.OffsetDateTime.now().format(java.time.format.DateTimeFormatter.ISO_DATE_TIME)}"

    """
    mkdir -p logs/elprep

    touch ${prefix}.elfasta
    touch logs/elprep/elprep-${timestamp}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(elprep 2>&1 | head -n2 | tail -n1 |sed 's/^.*version //;s/ compiled.*\$//')
    END_VERSIONS
    """
}
