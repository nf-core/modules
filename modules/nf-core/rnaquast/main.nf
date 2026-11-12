// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RNAQUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b8/b808671479318b6f8bb8d5b0c105ecccdea485b1f380e667366fa9363e15e61a/data'
        : 'community.wave.seqera.io/library/rnaquast:2.3.0--9fcd20e07559eb30'}"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('rnaquast'), eval("rnaQUAST.py -h | grep -i 'rnaQUAST.py v' | sed -E 's/.*v\\.([0-9.]+).*/\\1/'"), topic: versions, emit: versions_rnaquast

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    reference = reference ? "--reference ${reference}" : ''
    gtf = gtf ? "--gtf ${gtf}" : ''
    """
    # The gmap launcher exits 255 if any PATH entry is unreadable, e.g. /root/.pixi/bin when running as non-root
    export PATH=\$(echo "\$PATH" | tr ':' '\\n' | grep -v '^/root/' | paste -sd: -)

    rnaQUAST.py \\
        $args \\
        --threads $task.cpus \\
        --transcripts ${fasta} \\
        ${reference} \\
        ${gtf} \\
        -o ${prefix}
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/rnaQUAST.log
    """
}
