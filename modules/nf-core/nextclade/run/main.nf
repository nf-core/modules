process NEXTCLADE_RUN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/936786744b34cf016b948026a6b4e9489011424e15c28dfb2f7d03c31bb4afb5/data'
        : 'community.wave.seqera.io/library/nextclade:3.11.0--155203da8341cfe6'}"

    input:
    tuple val(meta), path(fasta)
    path dataset

    output:
    tuple val(meta), path("${prefix}.csv"), optional: true, emit: csv
    tuple val(meta), path("${prefix}.errors.csv"), optional: true, emit: csv_errors
    tuple val(meta), path("${prefix}.insertions.csv"), optional: true, emit: csv_insertions
    tuple val(meta), path("${prefix}.tsv"), optional: true, emit: tsv
    tuple val(meta), path("${prefix}.json"), optional: true, emit: json
    tuple val(meta), path("${prefix}.auspice.json"), optional: true, emit: json_auspice
    tuple val(meta), path("${prefix}.ndjson"), optional: true, emit: ndjson
    tuple val(meta), path("${prefix}.aligned.fasta"), optional: true, emit: fasta_aligned
    tuple val(meta), path("*_translation.*.fasta"), optional: true, emit: fasta_translation
    tuple val(meta), path("${prefix}.nwk"), optional: true, emit: nwk
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nextclade \\
        run \\
        ${args} \\
        --jobs ${task.cpus} \\
        --input-dataset ${dataset} \\
        --output-all ./ \\
        --output-basename ${prefix} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}
