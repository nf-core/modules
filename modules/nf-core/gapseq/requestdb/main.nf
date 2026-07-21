process GAPSEQ_REQUESTDB {
    tag "$taxon"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/933e301b11c1ec1699da6382e9e35b0e4e31edb80763eb2fa1b69ad7d6d1e5c7/data'
:         'community.wave.seqera.io/library/gapseq:2.1.0--c32b876ebb5e5f5b' }"

    input:
    val(taxon)

    output:
    path("gapseq_db"), emit: db
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Download the reference sequence database into a local directory.
    # gapseq resolves -D with readlink -f and sets userdir=true, so no
    # \$HOME/.gapseq/seq fallback is used.
    gapseq \\
        update-sequences \\
        -t $taxon \\
        -D gapseq_db
    """

    stub:
    """
    mkdir -p gapseq_db
    touch gapseq_db/seq_${taxon}.fa
    """
}
