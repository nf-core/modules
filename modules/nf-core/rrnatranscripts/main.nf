process RRNATRANSCRIPTS {
    tag '$rrna'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    path(gtf)

    output:
    path("*rrna_intervals.gtf") , emit: rrna_gtf, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "get_rrna_transcripts.py"

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"
    """
    touch ${prefix}_rrna_intervals.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
