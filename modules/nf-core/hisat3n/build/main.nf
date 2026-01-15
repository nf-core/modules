process HISAT3N_BUILD {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/01/019e3df43ceed54d16e3e97d86e448317fa7af5accb483763075d3ed2702737b/data' :
        'community.wave.seqera.io/library/hisat-3n:0.0.3--b4b80cb38c483147' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("hisat3n"), emit: index
    tuple val("${task.process}"), val('hisat-3n'), eval("hisat-3n --version 2>&1 | grep -o 'version [^ ]*' | cut -d ' ' -f 2"), topic: versions, emit: versions_hisat3n_build

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir hisat3n
    hisat-3n-build \\
        -p $task.cpus \\
        $args \\
        $fasta \\
        hisat3n/${fasta.baseName}
    """

    stub:
    """
    mkdir hisat3n
    """
}
