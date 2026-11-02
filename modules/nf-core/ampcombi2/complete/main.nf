// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AMPCOMBI2_COMPLETE {
    tag "ampcombi2"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/063288a19d4f267756c136b1ac8bd06ef89ecee62bf6fcc5b27732848f118b01/data'
        : 'community.wave.seqera.io/library/ampcombi:3.0.0--3b8b8b73f38675a2'}"

    input:
    path(summaries)

    output:
    path("Ampcombi_summary.tsv") , emit: tsv
    path("Ampcombi_complete.log"), emit: log, optional:true
    tuple val("${task.process}"), val('ampcombi'), eval("ampcombi --version | sed 's/ampcombi //'"), emit: versions_ampcombi, topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ampcombi complete \\
        --summaries_files '${summaries.collect{file_path -> "$file_path"}.join("' '")}' \\
        ${args}
    """

    stub:
    """
    touch Ampcombi_summary.tsv
    """
}
