process HUGGINGFACE_DOWNLOAD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/huggingface_hub:1.18.0--581812f40aa06159' :
        'community.wave.seqera.io/library/huggingface_hub:1.18.0--6bb77406fdf79177' }"

    input:
    tuple val(meta), val(hf_repo), val(hf_file)

    output:
    tuple val(meta), path(hf_file), emit: output
    tuple val("${task.process}"), val("huggingface_hub"), eval("hf --version 2>&1 | tail -n1 | awk '{print \$NF}'"), topic: versions, emit: versions_huggingface_hub


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
