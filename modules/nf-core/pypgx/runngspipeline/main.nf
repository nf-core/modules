process PYPGX_RUNNGSPIPELINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(coverage), path(control_stats), val(pgx_gene)
    tuple val(meta2), path(resource_bundle)

    output:
    tuple val(meta), path("*pypgx_output/results.zip"), emit: results
    tuple val(meta), path("*pypgx_output/cnv-calls.zip"), emit: cnv_calls, optional: true
    tuple val(meta), path("*pypgx_output/consolidated-variants.zip"), emit: consolidated_variants
    tuple val("${task.process}"), val('pypgx'), eval('pypgx -v 2>&1 | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+" | head -1'), emit: versions_pypgx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}_${pgx_gene}"
    def depth_coverage = coverage ? "--depth-of-coverage ${coverage}" : ""
    def control_statistics = control_stats ? "--control-statistics ${control_stats}" : ""

    """
    export MPLCONFIGDIR="/tmp/"
    export PYPGX_BUNDLE=${resource_bundle}/

    pypgx run-ngs-pipeline \\
        ${args} \\
        ${pgx_gene} \\
        ${prefix}_pypgx_output/ \\
        --variants ${vcf} \\
        ${depth_coverage} \\
        ${control_statistics}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${pgx_gene}"

    """
    mkdir ${prefix}_pypgx_output
    # zip program unavailable in container
    python -c 'import zipfile; zipfile.ZipFile("${prefix}_pypgx_output/results.zip", "w").close()'
    python -c 'import zipfile; zipfile.ZipFile("${prefix}_pypgx_output/cnv-calls.zip", "w").close()'
    python -c 'import zipfile; zipfile.ZipFile("${prefix}_pypgx_output/consolidated-variants.zip", "w").close()'
    """
}
