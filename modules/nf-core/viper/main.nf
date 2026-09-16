process VIPER {
    tag "$meta_exp.id vs $meta_net.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'papaemmelab/viper'

    input:
    tuple val(meta_net), path(network)
    tuple val(meta_exp), path(expression_matrix)

    output:
    path "${meta_exp.id}_vs_${meta_net.id}.viper_matrix.tsv", emit: viper_matrix
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'viper.R'

    stub:
    prefix = task.ext.prefix ?: "${meta_exp.id}_vs_${meta_net.id}"
    """
    touch ${prefix}.viper_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viper: \$(Rscript -e "cat(as.character(packageVersion('viper')))"
    END_VERSIONS
    """
}
