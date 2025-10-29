process KRAKEN2_BUILD {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data'
        : 'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71'}"

    input:
    tuple val(meta), path(library_added_files, stageAs: "database/library/added/")
    tuple val(meta2), path(seqid2taxid_map, stageAs: "database/seqid2taxid.map")
    tuple val(meta3), path(taxonomy_files, stageAs: "database/taxonomy/")
    val cleaning

    output:
    tuple val(meta), path("${prefix}"), emit: db
    tuple val(meta), path("${prefix}/*k2d"), path("${prefix}/*map"), path("${prefix}/library/added/*"), path("${prefix}/taxonomy/*"), optional: true, emit: db_separated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    run_clean = cleaning ? "kraken2-build --clean --db database/" : ""
    """
    kraken2-build \\
        --build \\
        ${args} \\
        --threads ${task.cpus} \\
        --db database/

    ${run_clean}

    mv database/ ${prefix}/

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
    mkdir -p "${prefix}"
    touch ${prefix}/{hash,opts,tax}.k2d

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
