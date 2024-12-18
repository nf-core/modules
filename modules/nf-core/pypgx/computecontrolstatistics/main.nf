process PYPGX_COMPUTECONTROLSTATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(control_gene)
    val(assembly_version)

    output:
    tuple val(meta), path("*.zip"), emit: control_stats
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = "${control_gene}"  ?: "VDR"
    def assembly = "${assembly_version}" ?: "GRCh38"

    """
    pypgx compute-control-statistics \\
        ${args} \\
        --assembly ${assembly} \\
        ${control} \\
        ${prefix}_${control}.zip \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def control = "${control_gene}"  ?: "VDR"
    """
    # zip program unavailable in container
    python -c 'import zipfile; zipfile.ZipFile("${prefix}_${control}.zip", "w").close()'
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
