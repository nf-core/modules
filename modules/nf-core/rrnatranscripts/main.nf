process RRNATRANSCRIPTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ca/ca12480e865de90051752075ee6b31c83559fa37775901622ca685b776912bb0/data'
        : 'community.wave.seqera.io/library/python:3.12.2--37a94ef16a29fce0'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*_rrna_intervals.gtf'), emit: rrna_gtf
    path "versions.yml"                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'get_rrna_transcripts.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rrna_intervals.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
