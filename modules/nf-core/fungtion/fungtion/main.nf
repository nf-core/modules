process FUNGTION_FUNGTION {
    tag "$meta.id"
    label 'process_high'

    conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66bda82170f41e26e1918642db86b02c04c0ecd7dd7d91937bcddcfe59f1d9cd/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4f0bfaae28ae79137790f705159bc756986d566f57622439712f8f4e6118eb1e/data') :
        (task.accelerator ? 'community.wave.seqera.io/library/python_r-base_r-caret_r-e1071_pruned:6e5900b868c2177c' : 'community.wave.seqera.io/library/python_r-base_r-caret_r-e1071_pruned:bb6702457fde82f3') }"

    input:
    tuple val(meta), path(fasta)
    path(pretrain)

    output:
    tuple val(meta), path("${prefix}/${prefix}.csv")        , emit: predictions
    tuple val(meta), path("${prefix}/${prefix}_analysis")   , emit: analysis,  optional: true
    tuple val(meta), path("${prefix}/${prefix}.html")       , emit: html,      optional: true
    tuple val(meta), path("${prefix}/${prefix}_assets")     , emit: assets,    optional: true
    tuple val(meta), path("${prefix}/${prefix}_temp_folder"), emit: temp,      optional: true
    tuple val(meta), path("${prefix}.log")                  , emit: log
    tuple val("${task.process}"), val('fungtion'), eval("fungtion --version 2>&1 | head -1"), topic: versions, emit: versions_fungtion
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('cuda'), eval('python -c "import torch; print(torch.version.cuda or \'no CUDA available\')"'), topic: versions, emit: versions_cuda

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def device = task.accelerator ? "cuda" : "cpu"
    """
    # Without this, the pipeline exit status can become tee's exit status,
    # which can mask a failing fungtion command.
    set -o pipefail

    fungtion \\
        $args \\
        --fasta ${fasta} \\
        --output-dir . \\
        --prefix ${prefix} \\
        --pretrain ${pretrain}/esm1b_t33_650M_UR50S.pt \\
        --device ${device} \\
        | tee ${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "$args"

    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.csv
    touch ${prefix}.log

    if [[ "$args" != *"--skip-visualization"* ]]; then mkdir -p ${prefix}/${prefix}_analysis; fi
    if [[ "$args" == *"--html-report"* ]]; then mkdir -p ${prefix}/${prefix}_assets && touch ${prefix}/${prefix}.html; fi
    if [[ "$args" == *"--keep-temp"* ]]; then mkdir -p ${prefix}/${prefix}_temp_folder; fi
    """
}
