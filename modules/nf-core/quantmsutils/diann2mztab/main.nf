process QUANTMSUTILS_DIANN2MZTAB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quantms-utils:0.0.23--pyh7e72e81_0' :
        'biocontainers/quantms-utils:0.0.23--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(report), path(report_pg), path(report_pr), path("version/versions.yml"), path(exp_design), path(ms_information), path(fasta)

    output:
    tuple val(meta), path("*msstats_in.csv"), emit: out_msstats
    tuple val(meta), path("*triqler_in.tsv"), emit: out_triqler
    tuple val(meta), path("*.mzTab"), optional: true, emit: out_mztab
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    quantmsutilsc diann2mztab \\
        --folder ./ \\
        --exp_design ${exp_design} \\
        --diann_version ./version/versions.yml \\
        ${args} \\
        2>&1 | tee convert_report.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: \$(pip show quantms-utils | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch test_sample_msstats_in.csv
    touch test_sample_triqler_in.tsv
    touch test_sample.mzTab
    touch convert_report.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quantms-utils: \$(pip show quantms-utils | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """
}
