process OPENMS_DECOYDATABASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84ad5437e080c0b3ec0ce315a54c0bb62950277a4496f4d1be911f82078dd2a1/data'
:         'community.wave.seqera.io/library/openms:3.6.0--88df48c0a426b7f7' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta"), emit: decoy_fasta
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}_decoy"

    """
    DecoyDatabase \\
        -in $fasta \\
        -out ${prefix}.fasta \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}_decoy"

    """
    touch ${prefix}.fasta
    """
}
