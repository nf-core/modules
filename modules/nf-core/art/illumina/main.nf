process ART_ILLUMINA {
    tag "${meta.id}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/art:2016.06.05--h589041f_9'
        : 'biocontainers/art:2016.06.05--h589041f_9'}"

    input:
    tuple val(meta), path(fasta)
    val sequencing_system
    val fold_coverage
    val read_length

    output:
    tuple val(meta), path("*.fq.gz"), emit: fastq
    tuple val(meta), path("*.aln"), emit: aln, optional: true
    tuple val(meta), path("*.sam"), emit: sam, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('art'), val("2016.06.05"), emit: versions_art_illumina, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    art_illumina \\
        -ss ${sequencing_system} \\
        -i ${fasta} \\
        -l ${read_length} \\
        -f ${fold_coverage} \\
        -o ${prefix} \\
        ${args}

    gzip \\
        --no-name \\
        ${args2} \\
        ${prefix}*.fq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fq.gz
    echo "" | gzip >  ${prefix}1.fq.gz
    echo "" | gzip >  ${prefix}2.fq.gz
    touch ${prefix}.aln
    touch ${prefix}1.aln
    touch ${prefix}2.aln
    touch ${prefix}.sam
    touch ${prefix}1.sam
    touch ${prefix}2.sam
    """
}
