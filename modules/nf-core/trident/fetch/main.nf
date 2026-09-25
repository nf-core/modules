process TRIDENT_FETCH {
    tag "${fetch_fn ?: ''} ${fetch_s ?: ''}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/poseidon-trident:2.2.2.1--hf7d7819_0'
        : 'quay.io/biocontainers/poseidon-trident:2.2.2.1--hf7d7819_0'}"

    input:
    tuple path(archive_dir, stageAs: "archive"), val(fetch_s), path(fetch_fn)

    output:
    // The local archive directory given as input, if any, is emitted as an output to allow for downstream processes to use it. Defaults to "output_archive" if no archive_dir is provided as input.
    path "${output_archive_dir}", type: 'dir', emit: local_archive, optional: true, includeInputs: true
    path "output_archive/*", type: 'dir', emit: downloaded_packages, optional: true
    tuple val("${task.process}"), val('trident'), eval('trident --version'), emit: versions_trident, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fetch_string = fetch_s ? "--fetchString ${fetch_s}" : ''
    def fetch_file = fetch_fn ? "--fetchFile ${fetch_fn}" : ''
    // fetch will always download to the first directory provided in `-d`, but check all provided dirs for already downloaded packages.
    // Handle multiple archive directories if provided
    def archives = archive_dir ? '-d ' + archive_dir.join(" -d ") : ''
    output_archive_dir = archive_dir ? archive_dir[0] : 'output_archive'
    """
    trident fetch\\
        -d output_archive/ \\
        ${archives} \\
        ${args} \\
        ${fetch_string} \\
        ${fetch_file}
    """

    stub:
    def args = task.ext.args ?: ''
    def fetch_string = fetch_s ? "--fetchString ${fetch_s}" : ''
    def fetch_file = fetch_fn ? "--fetchFile ${fetch_fn}" : ''
    def archives = archive_dir ? archive_dir.join(" -d ") : ''
    output_archive_dir = archive_dir ? archive_dir[0] : 'output_archive'
    """
    echo trident fetch ${archives} ${fetch_string} ${fetch_file} ${args}

    mkdir -p output_archive/dummy_package_dir
    touch output_archive/dummy_package_dir/POSEIDON.yml
    touch output_archive/dummy_package_dir/dummy_package.geno
    touch output_archive/dummy_package_dir/dummy_package.snp
    touch output_archive/dummy_package_dir/dummy_package.ind
    touch output_archive/dummy_package_dir/dummy_package.janno
    """
}
