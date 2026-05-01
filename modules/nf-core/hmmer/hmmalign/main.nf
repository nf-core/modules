process HMMER_HMMALIGN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07c4cbd91c4459dc86b13b5cd799cacba96b27d66c276485550d299c7a4c6f8a/data' :
        'community.wave.seqera.io/library/hmmer:3.4--cb5d2dd2e85974ca' }"

    input:
    tuple val(meta), path(fasta)
    path hmm

    output:
    tuple val(meta), path("*.sto.gz"), emit: sto
    tuple val("${task.process}"), val("hmmer"), eval("hmmalign -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER //'"), topic: versions, emit: versions_hmmer

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hmmalign \\
        $args \\
        $hmm \\
        $fasta | gzip -c > ${prefix}.sto.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.sto.gz
    """
}
