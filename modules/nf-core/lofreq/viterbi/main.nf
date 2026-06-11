process LOFREQ_VITERBI {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py310h47ef89e_10' :
        'quay.io/biocontainers/lofreq:2.1.5--py310h47ef89e_10' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('lofreq'), eval("lofreq version | awk 'NR==1{print \$2}'"), emit: versions_lofreq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq \\
        viterbi \\
        $args \\
        -ref $fasta \\
        $bam |
        samtools sort \\
            $args2 \\
            -T ${prefix} \\
            --threads $task.cpus \\
            -o ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
