process GFASTATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.9--h077b44d_1':
        'biocontainers/gfastats:1.3.9--h077b44d_1' }"

    input:
    tuple val(meta), path(assembly)
    val out_fmt
    val genome_size
    val target
    path agpfile
    path include_bed
    path exclude_bed
    path instructions

    output:
    tuple val(meta), path("*.assembly_summary"), emit: assembly_summary
    tuple val(meta), path("*.${out_fmt}.gz")   , emit: assembly
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def agp  = agpfile ? "--agp-to-path $agpfile" : ""
    def ibed = include_bed ? "--include-bed $include_bed" : ""
    def ebed = exclude_bed ? "--exclude-bed $exclude_bed" : ""
    def sak  = instructions ? "--swiss-army-knife $instructions" : ""
    """
    gfastats \\
        $args \\
        --threads $task.cpus \\
        $agp \\
        $ibed \\
        $ebed \\
        $sak \\
        --out-format ${prefix}.${out_fmt}.gz \\
        $assembly \\
        $genome_size \\
        $target \\
        > ${prefix}.assembly_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.${out_fmt}.gz
    touch ${prefix}.assembly_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """
}
