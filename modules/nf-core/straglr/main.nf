process STRAGLR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/98/985911fb3967ec3264e5447b609c521a57f217b3483464757b310c589057166f/data' :
        'community.wave.seqera.io/library/straglr:1.5.6--a02df76bc693eff5' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference), path(index)
    tuple val(meta3), path(loci)
    tuple val(meta4), path(exclude)
    tuple val(meta5), path(regions)
    

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('straglr'), eval("straglr.py --version"), emit: versions_straglr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def loci_arg    = loci    ? "--loci ${loci}"       : ''
    def exclude_arg = exclude ? "--exclude ${exclude}" : ''
    def regions_arg = regions ? "--regions ${regions}" : ''
    """
    straglr.py \\
        ${bam} \\
        ${reference} \\
        ${prefix} \\
        --nprocs ${task.cpus} \\
        ${loci_arg} \\
        ${exclude_arg} \\
        ${regions_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.bed 
    touch ${prefix}.vcf
    """
}
