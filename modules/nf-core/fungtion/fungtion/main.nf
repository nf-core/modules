process FUNGTION_FUNGTION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/401a8c91682bec5f7030ed167ab55986d9255b7ebde4b9b06ab6e5b1f447bd50/data':
        'community.wave.seqera.io/library/r-base_r-e1071_r-caret_r-optparse_pruned:c4d8b270ddf00cb9' }"

    input:
    tuple val(meta), path(fasta)
    path(pretrain)

    output:
    tuple val(meta), path("${prefix}/${prefix}.csv")             , emit: predictions
    tuple val(meta), path("${prefix}/${prefix}_analysis")        , emit: analysis,  optional: true
    tuple val(meta), path("${prefix}/${prefix}.html")            , emit: html,      optional: true
    tuple val(meta), path("${prefix}/${prefix}_assets")          , emit: assets,    optional: true
    tuple val(meta), path("${prefix}/${prefix}_temp_folder")     , emit: temp,      optional: true
    tuple val(meta), path("${prefix}.log")                       , emit: log
    tuple val("${task.process}"), val('fungtion'), eval("fungtion --version 2>&1 | head -1"), topic: versions, emit: versions_fungtion
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fungtion \\
        $args \\
        --fasta ${fasta} \\
        --output-dir . \\
        --prefix ${prefix} \\
        --pretrain ${pretrain}/esm1b_t33_650M_UR50S.pt \\
        1> ${prefix}.log
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
