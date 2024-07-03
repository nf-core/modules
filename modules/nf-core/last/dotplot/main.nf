process LAST_DOTPLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1542--h43eeafb_1' :
        'biocontainers/last:1542--h43eeafb_1' }"

    input:
    tuple val(meta), path(maf), path(annot_b)
    tuple val(meta2), path(annot_a)
    val(format)

    output:
    tuple val(meta), path("*.gif"), optional:true, emit: gif
    tuple val(meta), path("*.png"), optional:true, emit: png
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def annot_a_arg = annot_a ? "-a ${annot_a}" : ''
    def annot_b_arg = annot_b ? "-b ${annot_b}" : ''
    """
    last-dotplot \\
        $args \\
        $annot_a_arg \\
        $annot_b_arg \\
        $maf \\
        $prefix.$format

    # last-dotplot has no --version option so let's use lastal from the same suite
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version | sed 's/lastal //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch $prefix.$format

    # last-dotplot has no --version option so let's use lastal from the same suite
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version | sed 's/lastal //')
    END_VERSIONS
    """

}
