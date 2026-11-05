// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CHEWBBACA_CREATESCHEMA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/75d39fb0629de237b2e3f1c9faf7e093a3ca96a23bb94439dfe15a1f2bd6e552/data'
        : 'community.wave.seqera.io/library/chewbbaca_libmamba:de41ca94d8b05ad2'}"

    input:
    tuple val(meta), path(fasta, stageAs: "input_genomes/*")
    path prodigal_tf
    path cds

    output:
    tuple val(meta), path("results/$meta.id"), emit: schema
    path "results/cds_coordinates.tsv"       , emit: cds_coordinates
    path "results/invalid_cds.txt"           , emit: invalid_cds
    tuple val("${task.process}"), val("chewbbaca"), eval("chewie --version 2>&1 | sed 's/chewBBACA version: //'"), topic: versions, emit: versions_chewbbaca

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def schema = "--n ${prefix}"
    def prodigal_tf_opt = prodigal_tf ? "--ptf ${prodigal_tf}" : ""
    def cds_opt = cds ? "--cds ${cds}" : ""

    """
    find ./input_genomes/ -name "*.gz" | sed 's/.gz//' | xargs -I {} bash -c 'gzip -cdf {}.gz > {}'

    chewie \\
        CreateSchema \\
        -i input_genomes/ \\
        -o results \\
        $schema \\
        $args \\
        $prodigal_tf_opt \\
        $cds_opt \\
        --cpu $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p results/${prefix}/short/
    touch results/${prefix}/contigs-protein{1..5}.fasta
    touch results/${prefix}/.genes_list
    touch results/${prefix}/.schema_config
    touch results/${prefix}/short/contigs-protein{1..5}_short.fasta
    touch results/cds_coordinates.tsv
    touch results/invalid_cds.txt
    """
}
