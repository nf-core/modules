// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process NEXTCLADE_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a2/a2e78f92daf0513c2c72dd8101d5280ceeb1feb53e164ed1268fee6513449b0c/data' :
        'community.wave.seqera.io/library/nextclade:3.29.1--792e8bf7ff96aca1' }"

    input:
    tuple val(meta), path(fasta)
    path dataset

    output:
    tuple val(meta), path("${prefix}.csv")           , optional:true, emit: csv
    tuple val(meta), path("${prefix}.errors.csv")    , optional:true, emit: csv_errors
    tuple val(meta), path("${prefix}.insertions.csv"), optional:true, emit: csv_insertions
    tuple val(meta), path("${prefix}.tsv")           , optional:true, emit: tsv
    tuple val(meta), path("${prefix}.json")          , optional:true, emit: json
    tuple val(meta), path("${prefix}.auspice.json")  , optional:true, emit: json_auspice
    tuple val(meta), path("${prefix}.ndjson")        , optional:true, emit: ndjson
    tuple val(meta), path("${prefix}.aligned.fasta") , optional:true, emit: fasta_aligned
    tuple val(meta), path("*_translation.*.fasta")   , optional:true, emit: fasta_translation
    tuple val(meta), path("${prefix}.nwk")           , optional:true, emit: nwk
    tuple val("${task.process}"), val('nextclade'), eval("nextclade --version 2>&1 | sed 's/.*nextclade \\([^ ]*\\).*/\\1/'"), emit: versions_nextclade, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nextclade \\
        run \\
        $args \\
        --jobs $task.cpus \\
        --input-dataset $dataset \\
        --output-all ./ \\
        --output-basename ${prefix} \\
        $fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    touch ${prefix}.tsv
    touch ${prefix}.json
    touch ${prefix}.auspice.json
    touch ${prefix}.aligned.fasta
    touch ${prefix}.cds_translation.test.fasta
    touch ${prefix}.nwk
    """
}
