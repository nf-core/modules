process QUANTMSUTILS_MZMLSTATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quantms-utils:0.0.23--pyh7e72e81_0' :
        'biocontainers/quantms-utils:0.0.23--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(ms_file)

    output:
    tuple val(meta), path("*_ms_info.parquet"), emit: ms_statistics
    tuple val(meta), path("*_ms2_info.parquet"), emit: ms2_statistics, optional: true
    tuple val(meta), path("*_feature_info.parquet"), emit: feature_statistics, optional: true
    path "versions.yml", emit: versions
    tuple val(meta), path("*.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    quantmsutilsc mzmlstats --ms_path "${ms_file}" \\
        ${args} \\
        2>&1 | tee ${ms_file.baseName}_mzml_statistics.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: \$(pip show quantms-utils | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch "${meta.id}_ms_info.parquet"
    touch "${meta.id}_ms2_info.parquet"
    touch "${meta.id}_feature_info.parquet"
    touch "${ms_file.baseName}_mzml_statistics.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: 0.0.23
    END_VERSIONS
    """
}
