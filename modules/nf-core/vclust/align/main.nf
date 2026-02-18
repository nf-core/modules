process VCLUST_ALIGN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(filter)
    val save_alignment

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("*.ids.tsv")    , emit: ids
    tuple val(meta), path("*.aln.tsv")    , emit: aln, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def filter_argument = filter ? "--filter ${filter}" : ''
    save_alignment_option = save_alignment ? "--out-aln ${prefix}.aln.tsv" : ''
    """
    vclust \\
        align \\
        $args \\
        ${filter_argument} \\
        ${save_alignment_option} \\
        -t $task.cpus \\
        -i ${fasta} \\
        -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.ids.tsv
    touch ${prefix}.aln.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vclust: \$(vclust --version)
    END_VERSIONS
    """
}
