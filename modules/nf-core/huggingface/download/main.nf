process HUGGINGFACE_DOWNLOAD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ede89b6c560efe66ea1fc29f93947d4f9fe9f529de23929fb5cc55394c53bf5c/data'
        : 'community.wave.seqera.io/library/huggingface_hub:1.18.0--d5830d12561fd965' }"

    input:
    tuple val(meta), val(hf_repo), val(hf_file)

    output:
    tuple val(meta), path(output_path), emit: output
    tuple val("${task.process}"), val("huggingface_hub"), eval("hf --version 2>&1 | tail -n1 | awk '{print \$NF}'"), topic: versions, emit: versions_huggingface_hub

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_arg = hf_file ? "\"${hf_file}\"" : ''
    def local_dir = hf_file ? '.' : prefix
    output_path = hf_file ?: prefix
    """
    export HF_HOME="\$PWD/.hf_home"
    # Xet needs an explicit CA bundle in containers without a system trust store.
    export SSL_CERT_FILE="\${SSL_CERT_FILE:-\$(python -m certifi)}"

    hf download \\
        "${hf_repo}" \\
        ${file_arg} \\
        --local-dir "${local_dir}" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_path = hf_file ?: prefix
    def create_output = hf_file ? "touch \"${output_path}\"" : "mkdir -p \"${output_path}\""
    """
    mkdir -p "\$(dirname "${output_path}")"
    ${create_output}
    """
}
