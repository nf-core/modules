process CAALM_CAALM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/840046fbd06b709533c0f9443a9ab04663012feb074ebca60067c0adb76baa21/data':
        'community.wave.seqera.io/library/faiss-cpu_python_pip_caalm_torch:c3008a34cb7c94b7' }"

    input:
    tuple val(meta), path(fasta)
    tuple path(level0), path(level1), path(level2)

    output:
    tuple val(meta), path("${prefix}_predictions.tsv")      , emit: predictions
    tuple val(meta), path("${prefix}_probabilities.jsonl")  , emit: probabilities
    tuple val(meta), path("${prefix}_statistics.tsv")       , emit: statistics
    tuple val(meta), path("${prefix}_level0_embeddings.npy"), emit: embeddings_level0, optional: true
    tuple val(meta), path("${prefix}_level1_embeddings.npy"), emit: embeddings_level1, optional: true
    tuple val(meta), path("${prefix}_level2_embeddings.npy"), emit: embeddings_level2, optional: true
    tuple val(meta), path("${prefix}.log")                  , emit: log
    tuple val("${task.process}"), val('caalm'), eval("caalm --version 2>&1 | head -1"), topic: versions, emit: versions_caalm
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('torch'), eval("python -c 'import torch; print(torch.__version__)'"), topic: versions, emit: versions_torch
    tuple val("${task.process}"), val('faiss'), eval("python -c 'import faiss; print(faiss.__version__)'"), topic: versions, emit: versions_faiss

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    caalm \\
        $args \\
        --num-workers ${task.cpus} \\
        --level0-model ${level0} \\
        --level1-model ${level1} \\
        --level2-model ${level2}/model.pt \\
        --level2-faiss-dir ${level2}/faiss \\
        --level2-label-tsv-dir ${level2}/refdb \\
        --output-name ${prefix} \\
        -o . \\
        ${fasta} \\
        > ${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "$args"

    touch ${prefix}_predictions.tsv
    touch ${prefix}_probabilities.jsonl
    touch ${prefix}_statistics.tsv
    touch ${prefix}.log

    if [[ "$args" == *"--save-level0-embeddings"* ]]; then touch ${prefix}_level0_embeddings.npy; fi
    if [[ "$args" == *"--save-level1-embeddings"* ]]; then touch ${prefix}_level1_embeddings.npy; fi
    if [[ "$args" == *"--save-level2-embeddings"* ]]; then touch ${prefix}_level2_embeddings.npy; fi
    """
}
