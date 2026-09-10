process MUTSERVE_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c5f31fa5275ec069d6087fbd4589edce2613e8f02afcc5f6e00788af22e7355e/data'
        : 'community.wave.seqera.io/library/mutserve:2.0.3--9082a0e022e8ba1b'}"

    input:
    tuple val(meta), path(reads), path(index)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("*.txt"), emit: raw, optional: true
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true


    tuple val("${task.process}"), val('mutserve'), eval("mutserve --version 2>&1 | grep -o -E 'v[0-9.]+' | tr -d 'v'"), topic: versions, emit: versions_mutserve

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mutserve call \\
        ${args} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        --output ${prefix}.vcf.gz \\
        ${reads}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.txt
    touch ${prefix}.fasta
    """
}
