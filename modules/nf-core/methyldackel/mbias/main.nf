process METHYLDACKEL_MBIAS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7' :
        'quay.io/biocontainers/methyldackel:0.6.1--he4a0461_7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.mbias.txt"), emit: txt
    tuple val("${task.process}"), val('methyldackel'), eval("MethylDackel --version 2>&1 | cut -f1 -d' '"), emit: versions_methyldackel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MethylDackel mbias \\
        ${args} \\
        ${fasta} \\
        ${bam} \\
        ${prefix} \\
        --txt \\
        > ${prefix}.mbias.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mbias.txt
    """
}
