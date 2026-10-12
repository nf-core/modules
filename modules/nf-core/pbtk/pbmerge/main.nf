process PBTK_PBMERGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pbtk:3.1.1--h9ee0642_0'
        : 'quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0'}"

    input:
    tuple val(meta), path(bams, stageAs: "input/*")

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.pbi"), emit: pbi, optional: true
    tuple val("${task.process}"), val('pbtk'), eval("pbmerge --version | head -n1 | sed 's/pbmerge //' | sed -E 's/ .+//'"), emit: versions_pbmerge, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbmerge \\
        -o ${prefix}.bam \\
        -j ${task.cpus} \\
        ${args} \\
        ${bams}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.pbi
    """
}
