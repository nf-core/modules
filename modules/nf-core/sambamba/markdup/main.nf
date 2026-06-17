process SAMBAMBA_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d92b95f4b1fcff268d632d73b9adc861b0e2db41d4ac5ec1ae598f72f194b8fe/data':
        'community.wave.seqera.io/library/sambamba:1.0.1--f6f871dbcf29d001' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('sambamba'), eval("sambamba --version 2>&1 | grep -oPm1 'sambamba \\\\K[0-9.]+'"), topic: versions, emit: versions_sambamba


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sambamba \\
        markdup \\
        $args \\
        -t $task.cpus \\
        --tmpdir ./ \\
        $bam \\
        ${prefix}.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
