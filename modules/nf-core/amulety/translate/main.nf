process AMULETY_TRANSLATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/amulety_igblast_wget:7020cde3b45925f9':
        'community.wave.seqera.io/library/amulety_igblast_wget:e477bc17f7c35e7c' }"

    input:
    tuple val(meta), path(tsv)
    path(reference_igblast)

    output:
    tuple val(meta), path("*_translated.tsv"), emit: repertoire_translated
    tuple val("${task.process}"), val('amulety'), eval("amulety --help 2>&1 | grep -o 'version [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_amulety, topic: versions
    tuple val("${task.process}"), val('igblastn'), eval("igblastn -version | grep -o 'igblast[0-9\\. ]\\+' | grep -o '[0-9\\. ]\\+'"), emit: versions_igblastn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export IGDATA=${reference_igblast}
    amulety \\
    translate-igblast \\
    $args \\
    --input-file $tsv \\
    --output-dir . \\
    --reference-dir ${reference_igblast}

    mv *_translated.tsv ${prefix}_translated.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_translated.tsv
    """
}
