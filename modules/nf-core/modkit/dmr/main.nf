process MODKIT_DMR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ont-modkit:0.6.1--hcdda2d0_0'
        : 'quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0'}"

    input:
    tuple val(meta), path(bedmethyl_a), path(bedmethyl_a_tbi)
    tuple val(meta2), path(bedmethyl_b), path(bedmethyl_b_tbi)
    tuple val(meta3), path(regions_bed)
    tuple val(meta4), path(fasta)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('modkit'), eval("modkit --version | sed 's/modkit //'"), emit: versions_modkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-r ${regions_bed}" : ''
    """
    modkit \\
        dmr pair \\
        -a ${bedmethyl_a} \\
        -b ${bedmethyl_b} \\
        ${regions} \\
        --ref ${fasta} \\
        -o ${prefix}.bed \\
        --log-filepath ${prefix}.log \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.log
    """
}
