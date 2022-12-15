process EPANG_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::epa-ng=0.3.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epa-ng:0.3.8--h9a82719_1':
        'quay.io/biocontainers/epa-ng:0.3.8--h9a82719_1' }"

    input:
    tuple val(meta), path(refaln)
    path  fullaln

    output:
    tuple val(meta), path("query.fasta.gz"), emit: query
    path  "reference.fasta.gz"             , emit: reference
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    epa-ng \\
        $args \\
        --split $refaln $fullaln

    gzip query.fasta reference.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epang: \$(echo \$(epa-ng --version 2>&1) | sed 's/^EPA-ng v//')
    END_VERSIONS
    """
}
