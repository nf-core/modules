process PYPGX_PREPAREDEPTHOFCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(pgx_genes)
    val(assembly_version)

    output:
    tuple val(meta), path("*.zip"), emit: coverage
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes = "--genes ${pgx_genes.join(' ')}" ?: ''
    def assembly = assembly_version ?: "GRCh38"

    """
    pypgx prepare-depth-of-coverage \\
        ${args} \\
        ${genes} \\
        --assembly ${assembly} \\
        ${prefix}.zip \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c 'import zipfile; zipfile.ZipFile("${prefix}.zip", "w").close()'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
