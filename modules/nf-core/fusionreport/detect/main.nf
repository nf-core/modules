process FUSIONREPORT_DETECT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed32f46746a5b33a1b5d597cfe2f62d3b1cfa54638f57cadc5d5158d6a4959d2/data' :
        'community.wave.seqera.io/library/fusion-report_beautifulsoup4_click_colorlog_pruned:78488bd99166aa9a'}"

    input:
    tuple val(meta), path(arriba_fusions), path(starfusion_fusions), path(fusioncatcher_fusions)
    tuple val(meta2), path(fusionreport_ref)
    val(tools_cutoff)

    output:
    tuple val(meta), path("*fusionreport.tsv")           , emit: fusion_list
    tuple val(meta), path("*fusionreport_filtered.tsv")  , emit: fusion_list_filtered
    tuple val(meta), path("*index.html")                 , emit: report
    tuple val(meta), path("*_*.html")                    , emit: html                 , optional:true
    tuple val(meta), path("*.csv")                       , emit: csv                  , optional:true
    tuple val(meta), path("*.json")                      , emit: json                 , optional:true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def tools = arriba_fusions        ? "--arriba ${arriba_fusions} " : ''
    tools    += starfusion_fusions    ? "--starfusion ${starfusion_fusions} " : ''
    tools    += fusioncatcher_fusions ? "--fusioncatcher ${fusioncatcher_fusions} " : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusion_report run $meta.id . $fusionreport_ref $tools --allow-multiple-gene-symbols --tool-cutoff $tools_cutoff $args $args2

    mv fusion_list.tsv ${prefix}.fusionreport.tsv
    mv fusion_list_filtered.tsv ${prefix}.fusionreport_filtered.tsv
    mv index.html ${prefix}_fusionreport_index.html
    [ ! -f fusions.csv ] || mv fusions.csv ${prefix}.fusions.csv
    [ ! -f fusions.json ] || mv fusions.json ${prefix}.fusions.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
        fusion_report DB retrieval: \$(cat $fusionreport_ref/DB-timestamp.txt)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fusionreport_filtered.tsv
    touch ${prefix}.fusionreport.tsv
    touch ${prefix}_fusionreport_index.html
    touch AAA_BBB.html
    touch ${prefix}.fusions.csv
    touch ${prefix}.fusions.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
