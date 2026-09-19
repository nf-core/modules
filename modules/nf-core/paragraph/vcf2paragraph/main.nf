process PARAGRAPH_VCF2PARAGRAPH {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paragraph:2.3--h21f15d8_1':
        'quay.io/biocontainers/paragraph:2.3--h21f15d8_1' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.json.gz")  , emit: graph
    tuple val("${task.process}"), val('paragraph'), val('2.3'), emit: versions_paragraph, topic: versions
    tuple val("${task.process}"), val('bgzip'), eval("bgzip -h 2>&1 | tr '\n' ' ' | sed 's/^.*Version: //; s/ .*\$//'"), emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vcf2paragraph.py \\
        ${args} \\
        --reference-sequence ${fasta} \\
        ${vcf} \\
        ${prefix}.json

    bgzip --threads ${task.cpus} ${args2} ${prefix}.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.json.gz
    """
}
