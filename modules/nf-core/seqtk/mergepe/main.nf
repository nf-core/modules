process SEQTK_MERGEPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2dae00435c9dff9202076069bf0d252db26a3e8c708741f17b5fbd053886ff9a/data'
:         'community.wave.seqera.io/library/seqtk:1.4--3fa479d3961a5b8d' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: reads
    tuple val("${task.process}"), val('seqtk'), eval("seqtk 2>&1 | sed -n 's/^Version: //p'"), emit: versions_seqtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        ln -s ${reads} ${prefix}.fastq.gz
        """
    } else {
        """
        seqtk \\
            mergepe \\
            $args \\
            ${reads} \\
            | gzip -n >> ${prefix}.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
