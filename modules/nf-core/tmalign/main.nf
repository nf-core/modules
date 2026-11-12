process TMALIGN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/aec03e3d506eec5c0143081308c17360dacad0632b756252fe27d535250fbd38/data'
        : 'community.wave.seqera.io/library/tmalign:20240303--f2a25d77ead3fcd9'}"

    input:
    tuple val(meta), path(structures, stageAs: 'structures/*')
    tuple val(meta2), path(reference, stageAs: 'reference/*')

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('tmalign'), eval("TMalign -v | sed -n '/Version/{s/.*Version //;s/).*//;p}'"), topic: versions, emit: versions_tmalign

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls structures/ > structures.list

    # -dir1 aligns every listed structure against the reference in one run;
    # -outfmt 2 is TM-align's own tabular output. The trailing CPU-time
    # comment is dropped so the table is reproducible.
    TMalign \\
        -dir1 structures/ structures.list \\
        ${reference} \\
        -outfmt 2 \\
        ${args} \\
        | sed '/^#Total CPU time/d' \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
