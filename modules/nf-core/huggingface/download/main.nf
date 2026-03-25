process HF_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/huggingface_hub:1.6.0--c106a7f9664ca39b"

    input:
    tuple val(meta), val(hf_repo), val(hf_file)

    output:
    tuple val(meta), path(hf_file), emit: output
    tuple val("${task.process}"), val("huggingface_hub"), eval("hf --version 2>&1"), topic: versions, emit: versions_huggingface_hub

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export HF_HOME="${workflow.projectDir}/hf_cache"
    export HF_HUB_CACHE="${workflow.projectDir}/hf_cache"
    STARTDIR=\$PWD
    hf download ${hf_repo} ${hf_file} --local-dir \$STARTDIR
    """

    stub:
    """
    touch ${hf_file}
    """
}
