process SOURMASH_PAIRWISE {
    tag "$meta.id"
    label 'process_low'
    // Module structure based on nf-core/sourmash/compare

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/90c428a073e3c30e1f45fdc6e2ee39a9945682b50e71d347bd712ea59a9aa5be/data'
        : 'community.wave.seqera.io/library/sourmash_plugin_branchwater:0.9.14--3720b94abc654eac' }"

    input:
    tuple val(meta), path(signatures)
    path file_list

    output:
    tuple val(meta), path("*pairwise.csv"), emit: csv
    tuple val("${task.process}"), val("sourmash"), eval("sourmash --version 2>&1 | sed 's/^sourmash //'"), topic: versions, emit: versions_sourmash

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def combine_cmd = ''
    def pairwise_input = ''

    // input format handling
    if (file_list) {
        // file_list is a pathlist text file — use directly
        pairwise_input = file_list
    } else if (signatures instanceof Collection && signatures.size() > 1) {
        // multiple individual .sig files — combine into a zip first
        pairwise_input = "${prefix}_collection.zip"
        combine_cmd = "sourmash sig cat ${signatures.sort { it.toString() }.join(' ')} -o ${prefix}_collection.zip"
    } else {
        // single file (zip or sig) — use directly
        pairwise_input = signatures
    }

    """
    ${combine_cmd}
    sourmash scripts pairwise \\
        ${pairwise_input} \\
        -o ${prefix}_pairwise.csv \\
        --cores ${task.cpus} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}_pairwise.csv
    """
}
