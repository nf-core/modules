process POPSCLE_DSCPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popscle:0.1beta--h2c78cec_0' :
        'quay.io/biocontainers/popscle:0.1beta--h2c78cec_0' }"

    input:
    tuple val(meta), path(bam), path(vcf)

    output:
    tuple val(meta), path('*.cel.gz'), emit: cel
    tuple val(meta), path('*.plp.gz'), emit: plp
    tuple val(meta), path('*.var.gz'), emit: var
    tuple val(meta), path('*.umi.gz'), emit: umi
    tuple val("${task.process}"), val('popscle'), val('0.1'), emit: versions_popscle, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    popscle dsc-pileup \\
        --sam $bam \\
        --vcf $vcf \\
        --out $prefix \\
        $args \\
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.cel.gz
    echo "" | gzip > ${prefix}.var.gz
    echo "" | gzip > ${prefix}.plp.gz
    echo "" | gzip > ${prefix}.umi.gz
    """
}
