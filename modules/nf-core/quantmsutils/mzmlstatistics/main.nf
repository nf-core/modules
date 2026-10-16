process QUANTMSUTILS_MZMLSTATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7caaeee055202206dc635a6f3d154eb229162e924eda945fe76c88bda497aba/data'
:         'community.wave.seqera.io/library/quantms-utils_sdrf-pipelines_python_setuptools:a179872b9d0fa97d' }"

    input:
    tuple val(meta), path(ms_file)

    output:
    tuple val(meta), path("*_ms_info.parquet"), emit: ms_statistics
    tuple val(meta), path("*_ms2_info.parquet"), emit: ms2_statistics, optional: true
    tuple val(meta), path("*_feature_info.parquet"), emit: feature_statistics, optional: true
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('quantms-utils'), eval("quantmsutilsc --version 2>&1 | sed -n 's/quantmsutils //p'"), emit: versions_quantmsutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    quantmsutilsc mzmlstats \\
        --ms_path "${ms_file}" \\
        ${args} \\
        2>&1 | tee ${ms_file.baseName}_mzml_statistics.log
    """

    stub:
    """
    touch "${meta.id}_ms_info.parquet"
    touch "${meta.id}_ms2_info.parquet"
    touch "${meta.id}_feature_info.parquet"
    touch "${ms_file.baseName}_mzml_statistics.log"
    """
}
