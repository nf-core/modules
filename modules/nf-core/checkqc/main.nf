process CHECKQC {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ec3d6e7260c79ecd92ff53e66337a8f1db4ca8d0a3ba561c35f21fa9acdd6ba/data':
        'community.wave.seqera.io/library/sample-sheet_numpy_pandas_pip_pruned:0b9dc0869e46a949' }"

    input:
    tuple val(meta), path(run_dir)
    path(checkqc_config)

    output:
    tuple val(meta), path("*checkqc_report.json"), emit: report
    tuple val("${task.process}"), val('checkqc'), eval('checkqc --version | sed -e "s/checkqc, version //g"'), emit: versions_checkqc, topic: versions

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

    # Check if the output JSON file is empty
    if [[ ! -s checkqc_report.json ]] ; then
        echo "Error: Empty JSON files. Most likely due to missing files in run directory. See .command.log file for errors."
        exit 1
    fi
    """

    stub:
    """
    touch checkqc_report.json
    """
}
