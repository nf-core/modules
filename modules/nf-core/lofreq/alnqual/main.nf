process LOFREQ_ALNQUAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py310h47ef89e_10' :
        'quay.io/biocontainers/lofreq:2.1.5--py310h47ef89e_10' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('lofreq'), eval("lofreq version | awk 'NR==1{print \$2}'"), emit: versions_lofreq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    lofreq \\
        alnqual \\
        $args \\
        -b \\
        $bam \\
        $fasta > ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
