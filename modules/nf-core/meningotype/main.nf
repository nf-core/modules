process MENINGOTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::meningotype=0.8.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meningotype:0.8.5--pyhdfd78af_0' :
        'quay.io/biocontainers/meningotype:0.8.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    args += task.ext.custom_args ? ' ' + task.ext.custom_args : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meningotype \\
        $args \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meningotype: \$( echo \$(meningotype --version 2>&1) | sed 's/^.*meningotype v//' )
    END_VERSIONS
    """
}
