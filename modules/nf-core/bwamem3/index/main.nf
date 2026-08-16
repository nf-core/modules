process BWAMEM3_INDEX {
    tag "${meta.id}"
    // NOTE bwa-mem3 builds an FM-index with libsais; peak memory scales with the reference size.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt }

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/09/09c0e84ac22ba722bae6f75c5c450d72dfdc06e7224ad2c0bd042849632ca70e/data'
        : 'community.wave.seqera.io/library/bwa-mem3:0.8.0--a430e1e1df49a92e'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem3"), emit: index
    tuple val("${task.process}"), val('bwamem3'), eval("bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p'"), emit: versions_bwamem3, topic: versions

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
