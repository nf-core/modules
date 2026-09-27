process PBTK_PBINDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.1--h9ee0642_0':
        'quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pbi"), emit: pbi
    tuple val("${task.process}"), val('pbtk'), eval("pbindex --version | head -n1 | cut -d' ' -f2"), emit: versions_pbindex, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pbindex \\
        -j $task.cpus \\
        $bam
    """

    stub:
    """
    touch ${bam}.pbi
    """
}
