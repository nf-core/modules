process CHROMAP_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/809e5a0166357804ca10097100e96844100019a11a3aaebf14e2cceb2ee98c0a/data' :
        'community.wave.seqera.io/library/chromap:0.3.2--4ec4bca51cd82195' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.index"), emit: index
    tuple val("${task.process}"), val('chromap'), eval("chromap --version 2>&1"), topic: versions, emit: versions_chromap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    chromap \\
        -i \\
        $args \\
        -t $task.cpus \\
        -r $fasta \\
        -o ${prefix}.index
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.index
    """
}
