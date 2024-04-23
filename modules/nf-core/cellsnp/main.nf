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
    tuple val(meta), path("$prefix")     , emit: cellsnp_output
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
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp: \$(cellsnp-lite --v |& awk '{print \$2}')
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp: \$(cellsnp-lite --v |& awk '{print \$2}')
    END_VERSIONS
    """
}
