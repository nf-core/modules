process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::mosdepth=0.3.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.3--hdfd78af_1' :
        'quay.io/biocontainers/mosdepth:0.3.3--hdfd78af_1'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed
    val   window_size

    output:
    tuple val(meta), path('*.global.dist.txt')    , emit: global_txt
    tuple val(meta), path('*.region.dist.txt')    , emit: regions_txt , optional:true
    tuple val(meta), path('*.summary.txt')        , emit: summary_txt
    tuple val(meta), path('*.per-base.d4')        , emit: d4          , optional:true
    tuple val(meta), path('*.per-base.bed.gz')    , emit: per_base_bed, optional:true
    tuple val(meta), path('*.per-base.bed.gz.csi'), emit: per_base_csi, optional:true
    tuple val(meta), path('*.regions.bed.gz')     , emit: regions_bed , optional:true
    tuple val(meta), path('*.regions.bed.gz.csi') , emit: regions_csi , optional:true
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (window_size) {
        interval = "--by ${window_size}"
    } else if ( bed ) {
        interval = "--by ${bed}"
    } else {
        interval = ""
    }
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

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.global.dist.txt
    touch ${prefix}.region.dist.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.per-base.d4
    touch ${prefix}.per-base.bed.gz
    touch ${prefix}.per-base.bed.gz.csi
    touch ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}
