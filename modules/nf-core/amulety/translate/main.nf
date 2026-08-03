process AMULETY_TRANSLATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c39fc87288811f7806452ecbdb559b9e9bba71aebb82c60d60af939a73bdf614/data':
        'community.wave.seqera.io/library/igblast_curl_python_transformers_pruned:05685e2c81024d42' }"

    input:
    tuple val(meta), path(tsv)
    path(reference_igblast)

    output:
    tuple val(meta), path("*_translated.tsv"), emit: repertoire_translated
    tuple val("${task.process}"), val('amulety'), eval("amulety --help 2>&1 | grep -o 'version [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_amulety, topic: versions
    tuple val("${task.process}"), val('igblastn'), eval("igblastn -version | grep -o 'igblast[0-9\\. ]\\+' | grep -o '[0-9\\. ]\\+'"), emit: versions_igblastn, topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | grep -o 'Python [0-9\\.]\\+' | grep -o '[0-9\\.]\\+'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pytorch'), eval("python -c 'import torch; print(torch.__version__)'"), emit: versions_pytorch, topic: versions
    tuple val("${task.process}"), val('transformers'), eval("python -c 'import transformers; print(transformers.__version__)'"), emit: versions_transformers, topic: versions

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
