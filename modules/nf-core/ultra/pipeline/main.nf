process ULTRA_PIPELINE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.1--pyh7cba7a3_1':
        'quay.io/biocontainers/ultra_bioinformatics:0.1--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val("${task.process}"), val('ultra'), eval("uLTRA --version | sed 's/uLTRA //'"), emit: versions_ultra, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    uLTRA \\
        pipeline \\
        --t $task.cpus \\
        --prefix $prefix \\
        $args \\
        $fasta \\
        $gtf \\
        $reads \\
        ./
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sam
    """


}
