process BWAMEM3_INDEX {
    tag "$fasta"
    // NOTE bwa-mem3 builds an FM-index with libsais; peak memory scales with the reference size.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45a11b91903b7e31b9418ca5fe6852c82a2b1ee1b0772ad6623618cb0ba9da55/data'
        : 'community.wave.seqera.io/library/bwa-mem3:0.3.0--89c56b3ab74a5e5c'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem3"), emit: index
    tuple val("${task.process}"), val('bwamem3'), eval('bwa-mem3 version 2>/dev/null | head -n1 | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir bwamem3
    bwa-mem3 \\
        index \\
        $args \\
        -t ${task.cpus} \\
        -p bwamem3/${prefix} \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    """
    mkdir bwamem3
    touch bwamem3/${prefix}.0123
    touch bwamem3/${prefix}.amb
    touch bwamem3/${prefix}.ann
    touch bwamem3/${prefix}.bwt.2bit.64
    touch bwamem3/${prefix}.pac
    """
}
