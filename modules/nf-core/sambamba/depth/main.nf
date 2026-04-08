process SAMBAMBA_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d92b95f4b1fcff268d632d73b9adc861b0e2db41d4ac5ec1ae598f72f194b8fe/data':
        'community.wave.seqera.io/library/sambamba:1.0.1--f6f871dbcf29d001' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(bed)
    val(mode)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('sambamba'), eval("sambamba --version 2>&1 | grep -oPm1 'sambamba \\\\K[0-9.]+'"), topic: versions, emit: versions_sambamba
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!['region','window','base'].contains(mode)) {
        error "Mode needs to be one of: region, window, base"
    }
    def bed_arg = bed ? "--regions ${bed}" : ''

    """
    sambamba \\
        depth \\
        $mode \\
        $bed_arg \\
        $args \\
        -t $task.cpus \\
        -o ${prefix}.bed \\
        $bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
