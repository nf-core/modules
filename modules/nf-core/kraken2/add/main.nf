process KRAKEN2_ADD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data'
        : 'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71'}"

    input:
    tuple val(meta), path(fasta)
    path taxonomy_names, stageAs: 'taxonomy/names.dmp'
    path taxonomy_nodes, stageAs: 'taxonomy/nodes.dmp'
    path accession2taxid, stageAs: 'taxonomy/*'
    path seqid2taxid, stageAs: "seqid2taxid.map"

    output:
    tuple val(meta), path("${prefix}/library/added/*", includeInputs: true), emit: library_added_files
    tuple val(meta), path("${prefix}/seqid2taxid.map", includeInputs: true), optional: true, emit: seqid2taxid_map
    tuple val(meta), path("${prefix}/taxonomy/*", includeInputs: true), emit: taxonomy_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    inject_custom_seqid2taxid_map = seqid2taxid ? "cp ${seqid2taxid} ${prefix}/" : ""
    """
    mkdir -p ${prefix}
    mv "taxonomy" ${prefix}

    ${inject_custom_seqid2taxid_map}

    echo ${fasta} |\\
    tr -s " " "\\012" |\\
    xargs -I {} -n1 kraken2-build \\
        --add-to-library {} \\
        --db ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

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
    mkdir -p "${prefix}"/library/added "${prefix}"/taxonomy
    touch "${prefix}"/library/added/test.txt "${prefix}"/seqid2taxid.map "${prefix}"/taxonomy/{nodes,names}.dmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
