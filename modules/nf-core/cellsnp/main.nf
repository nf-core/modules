process CELLSNP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cellsnp-lite:1.2.3--h6141fd1_2' :
        'biocontainers/cellsnp-lite:1.2.3--h6141fd1_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(region_vcf), path(barcode), path(sample_list)

    output:
    tuple val(meta), path('*/cellSNP.base.vcf.gz')     , emit: base
    tuple val(meta), path('*/cellSNP.cells.vcf.gz')    , emit: cell, optional: true
    tuple val(meta), path('*/cellSNP.samples.tsv')             , emit: sample
    tuple val(meta), path('*/cellSNP.tag.AD.mtx')      , emit: allele_depth
    tuple val(meta), path('*/cellSNP.tag.DP.mtx')      , emit: depth_coverage
    tuple val(meta), path('*/cellSNP.tag.OTH.mtx')     , emit: depth_other
    path 'versions.yml'                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = barcode ? "-b $barcode" : ''
    def region_file = region_vcf ? "-R $region_vcf" : ''
    def sample = sample_list ? "-I $sample_list" : ''

    """
    cellsnp-lite -s $bam \\
        $barcode_file \\
        $region_file \\
        $sample \\
        -O $prefix \\
        --nproc $task.cpus \\
	--gzip \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp: \$(cellsnp-lite --v | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_file = barcode ? "-b $barcode" : ''
    def region_file = region_vcf ? "-R $region_vcf" : ''
    def sample = sample_list ? "-I $sample_list" : ''

    """
    mkdir $prefix
    touch $prefix/cellSNP.base.vcf.gz
    touch $prefix/cellSNP.samples.tsv
    touch $prefix/cellSNP.tag.AD.mtx
    touch $prefix/cellSNP.tag.DP.mtx
    touch $prefix/cellSNP.tag.OTH.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp: \$(cellsnp-lite --v | awk '{print \$2}')
    END_VERSIONS
    """
}
