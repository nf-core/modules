process GUNC_MERGECHECKM {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gunc=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gunc:1.0.5--pyhdfd78af_0' :
        'biocontainers/gunc:1.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gunc_file), path(checkm_file)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunc \\
        merge_checkm \\
        $args \\
        -g $gunc_file \\
        -c $checkm_file \\
        -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$( gunc --version )
    END_VERSIONS
    """
}
