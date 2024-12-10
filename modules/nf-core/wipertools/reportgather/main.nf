process WIPERTOOLS_REPORTGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.3--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("*_gathered.report") , emit: report_out
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    wipertools reportgather -r $reports -f ${prefix}_gathered.report ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools reportgather: \$(wipertools reportgather --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "GATHERED REPORT FROM MANY REPORTS" > ${prefix}_gathered.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools reportgather: \$(wipertools reportgather --version)
    END_VERSIONS
    """
}
