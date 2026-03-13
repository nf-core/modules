process STITCHR_THIMBLE {
    tag "$meta.id"
    label 'process_single'


    container "ghcr.io/qbic-pipelines/stitchr:1.3.1"


    input:
    tuple val(meta), path(samplesheet)
    tuple val(meta2), path(alt_codon_usage)
    tuple val(meta3), path(allele_preferences)
    val(species)


    output:
    tuple val(meta), path("${meta.id}.tsv"), emit: thimble
    tuple val("${task.process}"), val('stitchr'), eval("stitchr --version"), topic: versions, emit: versions_stitchr
    tuple val("${task.process}"), val('imgtgenedl'), val("0.6.1"), topic: versions, emit: versions_imgtgenedl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //optional files:
    def codon_usage = alt_codon_usage ? "-cu ${alt_codon_usage}" : ""
    def allele_preference = allele_preferences ? "-p ${allele_preferences}" : ""

    """
    thimble \
        $args \
        -in ${samplesheet} \
        -o ${prefix}.tsv \
        ${codon_usage} \
        ${allele_preference} \
        -s ${species}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.tsv
    """
}
