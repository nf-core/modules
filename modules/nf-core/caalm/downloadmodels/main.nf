process CAALM_DOWNLOADMODELS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81726ed58923a02e2e947be826509e0f34c3e762183e0bd5b6b109769b67ac84/data':
        'community.wave.seqera.io/library/huggingface_hub_python:91b414233adf3eb1' }"

    output:
    path("models"), emit: models
    tuple val("${task.process}"), val('huggingface_hub'), eval("python -c 'import huggingface_hub; print(huggingface_hub.__version__)'"), topic: versions, emit: versions_huggingface_hub

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Redirect all HuggingFace cache I/O to the task work directory —
    # \$HOME/.cache is read-only inside Singularity containers
    export HF_HOME=.
    export HF_HUB_DISABLE_XET=1
    python -c "from huggingface_hub import snapshot_download; snapshot_download('lczong/CAALM', local_dir='models')"
    """

    stub:
    """
    mkdir -p models/level0 models/level1 models/level2/faiss models/level2/refdb
    touch models/level2/model.pt
    """
}
