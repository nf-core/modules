process BISMARK_COVERAGE2CYTOSINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/38e61d14ccaed82f60c967132963eb467d0fa4bccb7a21404c49b4f377735f03/data' :
        'community.wave.seqera.io/library/bismark:0.25.1--1f50935de5d79c47' }"

    input:
    tuple val(meta), path(coverage_file)
    tuple val(meta2), path(fasta, stageAs: 'tmp/*') // This change mounts as directory containing the FASTA file to prevent nested symlinks
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.cov.gz")                      , emit: coverage,  optional: true
    tuple val(meta), path("*report.txt.gz")                , emit: report
    tuple val(meta), path("*cytosine_context_summary.txt") , emit: summary
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage2cytosine \\
        ${coverage_file} \\
        --genome ${index} \\
        --output ${prefix} \\
        --gzip \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.cov.gz
    echo | gzip > ${prefix}.report.txt.gz
    touch ${prefix}.cytosine_context_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
