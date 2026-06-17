process LOFREQ_INDELQUAL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4' :
        'quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('lofreq'), eval("lofreq version | sed -n '1s/^.* //p'"), emit: versions_lofreq, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq indelqual \\
        $args \\
        -f $fasta \\
        -o ${prefix}.bam \\
        $bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
