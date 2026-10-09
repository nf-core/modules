process KRAKEN2_ADD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ec5af88da27c32c52b667af2af67ce76ffc96715c1a40e0f4ab72dd87623adfa/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:b111fc3860b64b23' }"

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
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), topic: versions, emit: versions_kraken2

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
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"
    mkdir -p "${prefix}"/library/added "${prefix}"/taxonomy
    touch "${prefix}"/library/added/test.txt "${prefix}"/seqid2taxid.map "${prefix}"/taxonomy/{nodes,names}.dmp
    """
}
