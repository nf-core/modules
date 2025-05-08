process NGSBITS_BEDCOVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngs-bits:2025_03--py313h572c47f_0':
        'biocontainers/ngs-bits:2025_03--py313h572c47f_0' }"

    input:
    tuple val(meta), path(reads), path(reads_index), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def reference = fasta ? "-ref ${fasta}" : ""

    """
    BedCoverage \\
        $args \\
        -threads ${task.cpus} \\
        -bam ${reads} \\
        -in ${bed} \\
        -out ${prefix}.bed \\
        ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(BedCoverage --version 2>&1) | sed 's/BedCoverage //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngs-bits: \$(echo \$(BedCoverage --version 2>&1) | sed 's/BedCoverage //' )
    END_VERSIONS
    """
}
