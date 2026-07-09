process FUNGTION_DOWNLOADMODELS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/401a8c91682bec5f7030ed167ab55986d9255b7ebde4b9b06ab6e5b1f447bd50/data':
        'community.wave.seqera.io/library/r-base_r-e1071_r-caret_r-optparse_pruned:c4d8b270ddf00cb9' }"

    output:
    path("models")                                                                            , emit: models
    tuple val("${task.process}"), val('fungtion'), eval("fungtion --version 2>&1 | head -1"), topic: versions, emit: versions_fungtion
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    fungtion setup-models --model-dir models
    """

    stub:
    """
    mkdir -p models
    touch models/esm1b_t33_650M_UR50S.pt
    """
}
