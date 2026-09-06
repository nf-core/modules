process RMATS_PREP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rmats:4.3.0--py311hf2f0b74_5'
        : 'quay.io/biocontainers/rmats:4.3.0--py311hf2f0b74_5'}"

    input:
    tuple val(meta), path(genome_bam)
    tuple val(meta2), path(reference_gtf)
    val read_length

    output:
    tuple val(meta), path("*.rmats"), emit: rmats
    tuple val(meta), path("*read_outcomes_by_bam.txt"), emit: read_outcomes
    tuple val("${task.process}"), val('rmats'), eval('rmats.py --version | sed -e "s/v//g"'), emit: versions_rmats, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${genome_bam} > ${prefix}.prep.b1.txt

    rmats.py \\
        --task prep \\
        ${args} \\
        --nthread ${task.cpus} \\
        --b1 ${prefix}.prep.b1.txt \\
        --gtf ${reference_gtf} \\
        --readLength ${read_length} \\
        --tmp ${prefix}_rmats_tmp \\
        --od ${prefix}_rmats_prep

    cp ${prefix}_rmats_tmp/*.txt ${prefix}_read_outcomes_by_bam.txt
    cp ${prefix}_rmats_tmp/*.rmats ${prefix}.rmats
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.rmats
    touch ${prefix}_read_outcomes_by_bam.txt
    """
}
