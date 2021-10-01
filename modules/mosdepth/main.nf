// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::mosdepth=0.3.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mosdepth:0.3.2--h01d7912_0"
    } else {
        container "quay.io/biocontainers/mosdepth:0.3.2--h01d7912_0"
    }

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
    path  "versions.yml"                          , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def interval = window_size ? "--by ${window_size}" : "--by ${bed}"
    """
    mosdepth \\
        $interval \\
        $options.args \\
        $prefix \\
        $bam
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}
