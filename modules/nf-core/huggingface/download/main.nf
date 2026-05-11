process HUGGINGFACE_DOWNLOAD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff4a498ec5588fa68383dab4723c33f2f0fd1c737cd168a188d42b1bb229ecc5/data' :
        'community.wave.seqera.io/library/huggingface_hub:1.14.0--bcc2a8f142e78ee0' }"

    input:
    tuple val(meta), val(hf_repo), val(hf_file)

    output:
    tuple val(meta), path(hf_file), emit: output
    tuple val("${task.process}"), val("huggingface_hub"), eval("hf --version 2>&1 | tail -n1 | awk '{print \$NF}'"), topic: versions, emit: versions_huggingface_hub

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export HF_HOME="\$PWD/.hf_home"

    hf download \\
        ${hf_repo} \\
        ${hf_file} \\
        --local-dir \$PWD \\
        ${args}
    """

    stub:
    """
    touch ${hf_file}
    """
}
