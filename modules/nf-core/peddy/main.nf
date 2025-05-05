process PEDDY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/peddy:0.4.8--pyh5e36f6f_0' :
        'biocontainers/peddy:0.4.8--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple val(meta2), path(ped)
    tuple val(meta3), path(sites)

    output:
    tuple val(meta), path("${prefix}.vs.html")              , emit: vs_html
    tuple val(meta), path("${prefix}.html")                 , emit: html
    tuple val(meta), path("*.peddy.ped")                    , emit: ped
    tuple val(meta), path("*.het_check.png")                , optional: true, emit: het_check_png
    tuple val(meta), path("*.ped_check.png")                , optional: true, emit: ped_check_png
    tuple val(meta), path("*.sex_check.png")                , optional: true, emit: sex_check_png
    tuple val(meta), path("*.het_check.csv")                , optional: true, emit: het_check_csv
    tuple val(meta), path("*.ped_check.csv")                , optional: true, emit: ped_check_csv
    tuple val(meta), path("*.sex_check.csv")                , optional: true, emit: sex_check_csv
    tuple val(meta), path("*.ped_check.rel-difference.csv") , optional: true, emit: ped_check_rel_difference_csv
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sites_arg = sites ? "--sites $sites" : ''
    if (sites && args.contains('--sites')) error "Double definition of --sites (in sites channel and in ext.args)"
    """
    peddy \\
        $args \\
        --prefix $prefix \\
        --plot \\
        -p $task.cpus \\
        $vcf \\
        $sites_arg \\
        $ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peddy: \$( peddy --version 2>&1 | tail -1 | sed 's/peddy, version //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (sites && args.contains('--sites')) error "Double definition of --sites (in sites channel and in ext.args)"
    """
    touch ${prefix}.vs.html
    touch ${prefix}.html
    touch ${prefix}.peddy.ped
    touch ${prefix}.het_check.csv
    touch ${prefix}.ped_check.csv
    touch ${prefix}.sex_check.csv
    touch ${prefix}.het_check.png
    touch ${prefix}.ped_check.png
    touch ${prefix}.sex_check.png
    touch ${prefix}.ped_check.rel-difference.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peddy: \$( peddy --version 2>&1 | tail -1 | sed 's/peddy, version //' )
    END_VERSIONS
    """
}
