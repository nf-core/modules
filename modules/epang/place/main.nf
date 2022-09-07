process EPANG_PLACE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::epa-ng=0.3.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epa-ng:0.3.8--h9a82719_1':
        'quay.io/biocontainers/epa-ng:0.3.8--h9a82719_1' }"

    input:
    tuple val(meta), path(queryaln), path(referencealn), path(referencetree)

    output:
    tuple val(meta), path("*.epa_result.jplace"), emit: jplace
    path "*.epa_info.log"                       , emit: log
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    epa-ng \\
        $args \\
        --query $queryaln \\
        --ref-msa $referencealn \\
        --tree $referencetree \\
        --model $meta.model

    mv epa_result.jplace ${prefix}.epa_result.jplace
    mv epa_info.log ${prefix}.epa_info.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epang: \$(echo \$(epa-ng --version 2>&1) | sed 's/^EPA-ng v//')
    END_VERSIONS
    """
}
