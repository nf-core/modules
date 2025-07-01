process KRAKEN2_ADD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/29ed8f68315625eca61a3de9fcb7b8739fe8da23f5779eda3792b9d276aa3b8f/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:45764814c4bb5bf3' }"

    input:
    tuple val(meta), path(fasta)
    path taxonomy_names, stageAs: 'taxonomy/names.dmp'
    path taxonomy_nodes, stageAs: 'taxonomy/nodes.dmp'
    path accession2taxid, stageAs: 'taxonomy/*'
    path seqid2taxid, stageAs: "seqid2taxid.map"

    output:
    tuple val(meta), path("${prefix}"), emit: db
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
    mkdir "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
