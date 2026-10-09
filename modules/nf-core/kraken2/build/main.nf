process KRAKEN2_BUILD {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ec5af88da27c32c52b667af2af67ce76ffc96715c1a40e0f4ab72dd87623adfa/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:b111fc3860b64b23' }"

    input:
    tuple val(meta), path(library_added_files, stageAs: "kraken2-database/library/added/")
    tuple val(meta2), path(seqid2taxid_map, stageAs: "kraken2-database/seqid2taxid.map")
    tuple val(meta3), path(taxonomy_files, stageAs: "kraken2-database/taxonomy/")
    val cleaning

    output:
    tuple val(meta), path("kraken2-database"), emit: db
    tuple val(meta), path("kraken2-database/*k2d", includeInputs: true), path("kraken2-database/*map", includeInputs: true), path("kraken2-database/library/added/*", includeInputs: true), path("kraken2-database/taxonomy/*", includeInputs: true), optional: true, emit: db_separated
    tuple val(meta), path("kraken2-database/unmapped*.txt"), optional: true, emit: unmapped
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), topic: versions, emit: versions_kraken2

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
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p kraken2-database/
    touch kraken2-database/{hash,opts,tax}.k2d
    touch kraken2-database/unmapped.txt
    """
}
