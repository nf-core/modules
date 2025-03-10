process LAST_DOTPLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a35d17772276115874eda45f17b05407c9855ecec8c11d39cce49f0028f44cfb/data'
        : 'community.wave.seqera.io/library/last_open-fonts:f1048938520a62ad'}"

    input:
    tuple val(meta), path(maf), path(annot_b)
    tuple val(meta2), path(annot_a)
    val(format)
    val(filter)

    output:
    tuple val(meta), path("*.gif"), optional:true, emit: gif
    tuple val(meta), path("*.png"), optional:true, emit: png
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def annot_a_arg = annot_a ? "-a ${annot_a}" : ''
    def annot_b_arg = annot_b ? "-b ${annot_b}" : ''
    def input_command = filter ? "maf-linked ${args2}" : "zcat -f"
    """
    TTF=/home/ubuntu/conda_pkgs_dir/open-fonts-0.7.0-1/fonts/open-fonts/DejaVuSansMono-Regular.ttf
    [ -e "\$TTF" ] || TTF="/opt/conda/fonts/open-fonts/DejaVuSansMono-Regular.ttf"
    $input_command $maf |
    last-dotplot \\
        -f \$TTF \\
        $args \\
        $annot_a_arg \\
        $annot_b_arg \\
        - \\
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
