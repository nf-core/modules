process CHECKQC {
    label 'process_single'

    container "community.wave.seqera.io/library/python_pip_interop_checkqc:d76c912c8fadc561"

    input:
    path(run_dir)
    path(checkqc_config)

    output:
    path "*checkqc_report.json", emit: report
    path "versions.yml"        , emit: versions

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
        $run_dir > checkqc_report.json || test -s "checkqc_report.json"

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
