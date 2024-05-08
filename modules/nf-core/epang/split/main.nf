process EPANG_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epa-ng:0.3.8--h9a82719_1':
        'biocontainers/epa-ng:0.3.8--h9a82719_1' }"

    input:
    tuple val(meta), path(refaln), path(fullaln)

    output:
    tuple val(meta), path("*query.fasta.gz")    , emit: query
    tuple val(meta), path("*reference.fasta.gz"), emit: reference
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    epa-ng \\
        $args \\
        --split $refaln $fullaln

    gzip -c query.fasta > ${prefix}.query.fasta.gz; rm query.fasta
    gzip -c reference.fasta > ${prefix}.reference.fasta.gz; rm reference.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epang: \$(echo \$(epa-ng --version 2>&1) | sed 's/^EPA-ng v//')
    END_VERSIONS
    """
}
