process METHYLDACKEL_EXTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7' :
        'quay.io/biocontainers/methyldackel:0.6.1--he4a0461_7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bedGraph") , optional: true, emit: bedgraph
    tuple val(meta), path("*.methylKit"), optional: true, emit: methylkit
    tuple val("${task.process}"), val('methyldackel'), eval("MethylDackel --version 2>&1 | cut -f1 -d' '"), emit: versions_methyldackel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    MethylDackel extract \\
        ${args} \\
        ${fasta} \\
        ${bam}
    """

    stub:
    def args = task.ext.args ?: ''
    def out_extension = args.contains('--methylKit') ? 'methylKit' : 'bedGraph'
    """
    touch ${bam.baseName}_CpG.${out_extension}
    """
}
