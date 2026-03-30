process CAALM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/840046fbd06b709533c0f9443a9ab04663012feb074ebca60067c0adb76baa21/data':
        'community.wave.seqera.io/library/faiss-cpu_python_pip_caalm_torch:c3008a34cb7c94b7' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_predictions.tsv")      , emit: predictions
    tuple val(meta), path("${prefix}_probabilities.jsonl")  , emit: probabilities
    tuple val(meta), path("${prefix}_statistics.tsv")       , emit: statistics
    tuple val(meta), path("${prefix}_level0_embeddings.npy"), emit: embeddings_level0, optional: true
    tuple val(meta), path("${prefix}_level1_embeddings.npy"), emit: embeddings_level1, optional: true
    tuple val(meta), path("${prefix}_level2_embeddings.npy"), emit: embeddings_level2, optional: true
    tuple val("${task.process}"), val('caalm'), eval("caalm --version 2>&1 | head -1"), topic: versions, emit: versions_caalm

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    caalm \\
        $args \\
        --output-name ${prefix} \\
        -o . \\
        ${fasta}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "$args"

    touch ${prefix}_predictions.tsv
    touch ${prefix}_probabilities.jsonl
    touch ${prefix}_statistics.tsv

    if [[ "$args" == *"--save-level0-embeddings"* ]]; then touch ${prefix}_level0_embeddings.npy; fi
    if [[ "$args" == *"--save-level1-embeddings"* ]]; then touch ${prefix}_level1_embeddings.npy; fi
    if [[ "$args" == *"--save-level2-embeddings"* ]]; then touch ${prefix}_level2_embeddings.npy; fi
    """
}
