process CHECKQC {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/python_numpy_pip_checkqc_interop:b5301d9801b8e66b"

    input:
    tuple val(meta), path(run_dir)
    path(checkqc_config)

    output:
    tuple val(meta), path("*checkqc_report.json"), emit: report
    tuple val(meta), path("*checkqc_log.txt")    , emit: log
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CheckQC module does not support Conda yet. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def config = checkqc_config ? "--config $checkqc_config" : ''

    """
    checkqc \
        $args \
        $config \
        --json \
        $run_dir > checkqc_report.json || true

    cp .command.log checkqc_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkqc: \$( checkqc --version | sed -e "s/checkqc, version //g" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch checkqc_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkqc: \$( checkqc --version | sed -e "s/checkqc, version //g" )
    END_VERSIONS
    """
}
