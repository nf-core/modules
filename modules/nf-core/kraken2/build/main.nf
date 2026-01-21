process KRAKEN2_BUILD {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data'
        : 'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71'}"

    input:
    tuple val(meta), path(library_added_files, stageAs: "kraken2-database/library/added/")
    tuple val(meta2), path(seqid2taxid_map, stageAs: "kraken2-database/seqid2taxid.map")
    tuple val(meta3), path(taxonomy_files, stageAs: "kraken2-database/taxonomy/")
    val cleaning

    output:
    tuple val(meta), path("kraken2-database"), emit: db
    tuple val(meta), path("kraken2-database/*k2d", includeInputs: true), path("kraken2-database/*map", includeInputs: true), path("kraken2-database/library/added/*", includeInputs: true), path("kraken2-database/taxonomy/*", includeInputs: true), optional: true, emit: db_separated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    run_clean = cleaning ? "kraken2-build --clean --db kraken2-database/" : ""
    """
    kraken2-build \\
        --build \\
        ${args} \\
        --threads ${task.cpus} \\
        --db kraken2-database/

    ${run_clean}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p kraken2-database/
    touch kraken2-database/{hash,opts,tax}.k2d

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
