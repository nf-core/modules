process REGTOOLS_JUNCTIONSEXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/143bc2dfa40320fbe1e52329953c8508780591835223f4ca492d3206598604a8/data' :
        'community.wave.seqera.io/library/regtools:1.0.0--461ddf16709a70cf' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(strand_specificity)

    output:
    tuple val(meta), path("*.junc"), emit: junc
    tuple val("${task.process}"), val('regtools'), eval("regtools --version 2>&1 | sed -n 's/Version:\t//p'"), topic: versions, emit: versions_repeatmasker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strand = strand_specificity ?: 'XS'
    """
    regtools junctions extract \\
        $args \\
        -s ${strand} \\
        -o ${prefix}.junc \\
        $bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junc
    """
}
