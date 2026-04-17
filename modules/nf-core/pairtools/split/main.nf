process PAIRTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.1.3--py39h7a39fba_0' :
        'biocontainers/pairtools:1.1.3--py39h7a39fba_0' }"

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path("*.split.pairs.gz"), emit: pairs
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val("${task.process}"), val('pairtools'), eval("pairtools --version | sed 's/.*pairtools.*version //'") , emit: versions_pairtools, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=tmp
    echo "" | gzip > ${prefix}.split.pairs.gz
    """
}
