process STAPHSCAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staphscan:0.3.1--pyhdfd78af_0' :
        'quay.io/biocontainers/staphscan:0.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('staphscan'), eval("staphscan --version | sed 's/staphscan //;'"), emit: versions_staphscan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    staphscan \\
        $args \\
        -o ${prefix}_outdir \\
        --report ${prefix}_summary.tsv \\
        -i $fastas
    mv ${prefix}_outdir/*.tsv . || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.tsv
    """
}
