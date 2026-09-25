process NARFMAP_HASHTABLE {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e09c4813f50c84494bc9e482c8e5ac52a0c3d5567945b7056c8c36d54aff894/data':
        'community.wave.seqera.io/library/narfmap_pigz_samtools:e3bfa7f4d4cfb1bb' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("narfmap")    , emit: hashmap
    tuple val("${task.process}"), val('narfmap'), eval("dragen-os --version 2>&1"), topic: versions, emit: versions_narfmap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir narfmap
    dragen-os \\
        --build-hash-table true \\
        --ht-reference ${fasta} \\
        --output-directory narfmap \\
        $args \\
        --ht-num-threads ${task.cpus}
    """

    stub:
    """
    mkdir narfmap
    """
}
