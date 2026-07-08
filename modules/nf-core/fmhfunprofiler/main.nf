process FMHFUNPROFILER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fmh-funprofiler:1.1.1--pyh106432d_0':
        'biocontainers/fmh-funprofiler:1.1.1--pyh106432d_0' }"

    input:
    tuple val(meta), path(metagenome)
    path(ko_sketch)
    val(ksize)
    val(scaled)

    output:
    tuple val(meta), path("*.csv") , emit: csv
    tuple val(meta), path("*.biom"), emit: biom, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    funcprofiler \\
        $args \\
        $metagenome \\
        $ko_sketch \\
        $ksize \\
        $scaled \\
        ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fmh-funprofiler: \$(funcprofiler --version 2>&1 | sed 's/funcprofiler //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    touch ${prefix}.biom

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fmh-funprofiler: \$(funcprofiler --version 2>&1 | sed 's/funcprofiler //')
    END_VERSIONS
    """
}