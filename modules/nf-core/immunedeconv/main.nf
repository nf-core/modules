process IMMUNEDECONV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/22cb85f1b69ceff45b83e0fdb7b96d9ae29c8aafeaa0707d64cc4628982977ab/data' :
    'community.wave.seqera.io/library/r-immunedeconv:2.1.2--e1bb1ea1cf505cb3' }"

    input:
    tuple val(meta), path(input_file), val(method), val(function)
    val gene_symbol_col

    output:
    tuple val(meta), path("*.deconvolution_results.tsv"), emit: deconv_table
    tuple val(meta), path("*.png"),                       emit: deconv_plots, optional: true
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'immunedeconv.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.deconvolution_results.tsv
    touch ${prefix}.plot1_stacked_bar_chart.png
    touch ${prefix}.plot2_points_with_facets.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-immunedeconv: \$(Rscript -e "cat(as.character(packageVersion('immunedeconv')))")
    END_VERSIONS
    """
}
