process RIBOTISH_QUALITY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f77d9ee4f52ca5508563a49e52a80879d5b22e6a71824723f5865cad8c515d86/data'
        : 'community.wave.seqera.io/library/ribotish:0.2.8--46a4cf0c323417eb'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.txt")    , emit: distribution
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.para.py"), emit: offset
    tuple val("${task.process}"), val('ribotish'), eval("ribotish --version | sed 's/ribotish //'"), topic: versions, emit: versions_ribotish

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ribotish quality \\
        -b $bam \\
        -g $gtf \\
        -o ${prefix}_qual.txt \\
        -f ${prefix}_qual.pdf \\
        -r ${prefix}.para.py \\
        -p $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_qual.txt ${prefix}_qual.pdf ${prefix}.para.py
    """
}
