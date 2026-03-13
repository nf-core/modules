process LAST_DOTPLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386dbe9821e6f2eebd5bcc67b3c82e2b730135c51bb565bae6511b24aba17e56/data'
        : 'community.wave.seqera.io/library/last_open-fonts:b8d1af8fd12256e2'}"

    input:
    tuple val(meta), path(maf), path(annot_b)
    tuple val(meta2), path(annot_a)
    val(format)
    val(filter)

    output:
    tuple val(meta), path("*.{gif,png}"), emit: plot
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions

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
    TTF=/home/runner/conda_pkgs_dir/open-fonts-0.7.0-1/fonts/open-fonts/DejaVuSansMono-Regular.ttf
    [ -e "\$TTF" ] || TTF="/opt/conda/fonts/open-fonts/DejaVuSansMono-Regular.ttf"
    $input_command $maf |
    last-dotplot \\
        -f \$TTF \\
        $args \\
        $annot_a_arg \\
        $annot_b_arg \\
        - \\
        $prefix.$format
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch $prefix.$format
    """

}
