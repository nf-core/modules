process BISMARK_COVERAGE2CYTOSINE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe2a9d58209a38df5c99615a41d9ed6e8e546380d04c176e076e107010819a72/data' :
        'community.wave.seqera.io/library/bismark:0.25.0--95ba99b483e2eaf9' }"

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cov.gz
    touch ${prefix}.report.txt.gz
    touch ${prefix}.cytosine_context_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
