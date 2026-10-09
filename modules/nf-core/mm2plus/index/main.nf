process MM2PLUS_INDEX {
    label 'process_low'

    // Note: the versions here need to match the versions used in mm2plus/align
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/72/7224ad076c5513325c2ef76f8788249af1d791dacaf9ee378d5b7512019e3a32/data' :
        'community.wave.seqera.io/library/mm2plus_samtools:fe581c94b0a4dc10' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    tuple val("${task.process}"), val("mm2plus"), eval("mm2plus --version"), topic: versions, emit: versions_mm2plus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mm2plus \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta
    """

    stub:
    """
    touch ${fasta.baseName}.mmi
    """
}
