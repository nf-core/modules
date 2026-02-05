process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14e1d96665f934a98e569fc5a6fa237f98d3753eee2b6f60d0aea8ff9d44f406/data' :
        'community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(variant_catalog)

    output:
    tuple val(meta), path("*.vcf.gz")        , emit: vcf
    tuple val(meta), path("*.json.gz")       , emit: json
    tuple val(meta), path("*_realigned.bam") , emit: bam
    tuple val("${task.process}"), val('expansionhunter'), eval("ExpansionHunter --version | head -1 | sed -n 's/^.*ExpansionHunter v//; s/]//p'"), topic: versions, emit: versions_expansionhunter
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ExpansionHunter \\
        ${args} \\
        --reads ${bam} \\
        --output-prefix ${prefix} \\
        --reference ${fasta} \\
        --variant-catalog ${variant_catalog}

    bgzip --threads ${task.cpus} ${args2} ${prefix}.vcf
    bgzip --threads ${task.cpus} ${args2} ${prefix}.json

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.json.gz
    touch ${prefix}_realigned.bam

    """
}
