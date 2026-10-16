process QUANTMSUTILS_DIANN2MZTAB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7caaeee055202206dc635a6f3d154eb229162e924eda945fe76c88bda497aba/data'
:         'community.wave.seqera.io/library/quantms-utils_sdrf-pipelines_python_setuptools:a179872b9d0fa97d' }"

    input:
    tuple val(meta), path(report), path(report_pg), path(report_pr), val(diann_version), path(exp_design), path(ms_information), path(fasta)

    output:
    tuple val(meta), path("*msstats_in.csv"), emit: out_msstats
    tuple val(meta), path("*triqler_in.tsv"), emit: out_triqler
    tuple val(meta), path("*.mzTab"), optional: true, emit: out_mztab
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('quantms-utils'), eval("quantmsutilsc --version 2>&1 | sed -n 's/quantmsutils //p'"), emit: versions_quantmsutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir version
    cat <<-END_VERSIONS > ./version/versions.yml
    "${task.process}":
        DIA-NN: ${diann_version}
    END_VERSIONS

    quantmsutilsc diann2mztab \\
        --folder ./ \\
        --exp_design ${exp_design} \\
        --diann_version ./version/versions.yml \\
        ${args} \\
        2>&1 | tee convert_report.log
    """

    stub:
    """
    touch test_sample_msstats_in.csv
    touch test_sample_triqler_in.tsv
    touch test_sample.mzTab
    touch convert_report.log
    """
}
