process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::mosdepth=0.3.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.2--h01d7912_0' :
        'quay.io/biocontainers/mosdepth:0.3.2--h01d7912_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed
    val   window_size

    output:
    tuple val(meta), path('*.global.dist.txt')    , emit: global_txt
    tuple val(meta), path('*.region.dist.txt')    , emit: regions_txt
    tuple val(meta), path('*.summary.txt')        , emit: summary_txt
    tuple val(meta), path('*.per-base.bed.gz')    , emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi'), emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')     , emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi') , emit: regions_csi
    path  "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def interval = window_size ? "--by ${window_size}" : "--by ${bed}"
    """
    mosdepth \\
        $interval \\
        $args \\
        $prefix \\
        $bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}
