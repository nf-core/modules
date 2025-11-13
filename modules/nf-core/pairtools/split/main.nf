process PAIRTOOLS_SPLIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pairtools:1.0.2--py39h2a9f597_0'
        : 'biocontainers/pairtools:1.0.2--py39h2a9f597_0'}"

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path("*.split.pairs.gz"), emit: pairs
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=tmp
    pairtools split \
        --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
        --output-pairs ${prefix}.split.pairs.gz \
        ${args} \
        ${pairs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=tmp
    echo "" | gzip > ${prefix}.split.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | sed 's/pairtools, version //')
    END_VERSIONS
    """
}
