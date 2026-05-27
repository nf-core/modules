process GAPSEQ_FILL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gapseq:2.0.1--hdfd78af_0' :
        'quay.io/biocontainers/gapseq:2.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(draft), path(medium)

    output:
    tuple val(meta), path("*-filled.RDS")  , emit: filled
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def medium_arg = medium ? "-n $medium" : "-n \$gapseq_root/dat/media/ALLmed.csv"
    """
    gapseq_bin=\$(readlink -f \$(which gapseq))
    gapseq_root=\$(dirname "\$gapseq_bin")

    draft_model=\$(ls *-draft.RDS)

    gapseq \\
        fill \\
        -m \$draft_model \\
        $medium_arg \\
        $args

    # Rename output file
    filled_model=\$(ls *.RDS | grep -Ev '(-draft)\\.RDS' | head -n1)
    mv "\$filled_model" ${prefix}-filled.RDS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-filled.RDS
    """
}
