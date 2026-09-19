process GAPSEQ_FILL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/933e301b11c1ec1699da6382e9e35b0e4e31edb80763eb2fa1b69ad7d6d1e5c7/data'
:         'community.wave.seqera.io/library/gapseq:2.1.0--c32b876ebb5e5f5b' }"

    input:
    tuple val(meta), path(draft), path(medium)

    output:
    tuple val(meta), path("*-filled.RDS")  , emit: filled
    tuple val(meta), path("*-filled.xml")  , emit: xml
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

    # Rename output files
    filled_model=\$(ls *.RDS | grep -Ev '(-draft)\\.RDS' | head -n1)
    mv "\$filled_model" ${prefix}-filled.RDS

    filled_xml=\$(ls *.xml | grep -Ev '(-draft)\\.xml' | head -n1)
    mv "\$filled_xml" ${prefix}-filled.xml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-filled.RDS
    touch ${prefix}-filled.xml
    """
}
