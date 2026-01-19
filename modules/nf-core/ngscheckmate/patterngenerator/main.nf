process NGSCHECKMATE_PATTERNGENERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a87c8e024fc7d44064c5c304d3b3bd668a88579b9e069d40b74bcc2458d9dc91/data':
        'community.wave.seqera.io/library/bcftools_ngscheckmate:5d4fe3e82ae99a2b' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bowtie_index)

    output:
    tuple val(meta), path("*.pt"), emit: pt
    tuple val("${task.process}"), val('ngscheckmate'), eval("ncm.py --help | sed '7!d;s/.* v//g'"), topic: versions, emit: versions_ngscheckmate

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "makesnvpattern.pl generates a fasta file with the same name as the input fasta, use \"task.ext.prefix\" to disambiguate!"

    """
    INDEX=\$(find -L ./ -name "*.3.ebwt" | sed 's/\\.3.ebwt\$//')

    makesnvpattern.pl ${bed} ${fasta} \$INDEX . ${prefix}
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "makesnvpattern.pl generates a fasta file with the same name as the input fasta, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.pt
    """
}
