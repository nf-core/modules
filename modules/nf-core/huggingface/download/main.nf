process HUGGINGFACE_DOWNLOAD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b2a780393f48a8b97b0510276123e5f2ff11eddf0a3382471f29c68320484e3/data' :
        'community.wave.seqera.io/library/huggingface_hub:1.6.0--c106a7f9664ca39b' }"

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
