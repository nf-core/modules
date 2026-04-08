process AMULETY_TRANSLATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/eef2fbc7c8d1ba71b3890a83b520c3eefa136ec4de5a8e6a97db828ae354d7ab/data':
        'community.wave.seqera.io/library/igblast_curl_python_wget_pruned:07dda71433b05ed5' }"

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
    --nproc ${task.cpus} \\
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
