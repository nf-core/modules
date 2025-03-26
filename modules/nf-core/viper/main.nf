process VIPER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'papaemmelab/viper'

    input:
    tuple val(meta), path(network)
    path expression_matrix

    output:
    path "viper_results.tsv", emit: viper_results
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'viper.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch viper_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viper: \$(Rscript -e "cat(as.character(packageVersion('viper')))"
    END_VERSIONS
    """
}
