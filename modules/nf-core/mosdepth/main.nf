process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00d32b53160c26794959da7303ee6e2107afd4d292060c9f287b0af1fddbd847/data' :
        'community.wave.seqera.io/library/mosdepth_htslib:0f58993cb6d93294'}"

    input:
    tuple val(meta),  path(bam), path(bai), path(bed)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    tuple val("${task.process}"), val('mosdepth'), eval("mosdepth --version | sed 's/mosdepth //g'"), topic: versions, emit: versions_mosdepth

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && (args.contains("--by") || args.contains("-b "))) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (args.contains("--thresholds") && !(bed || args.contains("--by") || args.contains("-b "))) {
        error "'--thresholds' can only be specified in conjunction with '--by' or an input bed file"
    }

    """
    mosdepth \\
        --threads $task.cpus \\
        $interval \\
        $reference \\
        $args \\
        $prefix \\
        $bam
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (bed && (args.contains("--by") || args.contains("-b "))) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (args.contains("--thresholds") && !(bed || args.contains("--by") || args.contains("-b "))) {
        error "'--thresholds' can only be specified in conjunction with '--by' or an input bed file"
    }
    """
    touch ${prefix}.global.dist.txt
    touch ${prefix}.region.dist.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.per-base.d4
    echo "" | gzip > ${prefix}.per-base.bed.gz
    touch ${prefix}.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.quantized.bed.gz
    touch ${prefix}.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.thresholds.bed.gz
    touch ${prefix}.thresholds.bed.gz.csi
    """
}
