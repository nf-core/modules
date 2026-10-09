process POPSCLE_DEMUXLET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popscle:0.1beta--h2c78cec_0' :
        'quay.io/biocontainers/popscle:0.1beta--h2c78cec_0' }"

    input:
    tuple val(meta), val(plp_prefix), path(bam), path(donor_genotype)

    output:
    tuple val(meta), path('*.best'), emit: demuxlet_result
    tuple val("${task.process}"), val('popscle'), val('0.1'), emit: versions_popscle, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = plp_prefix ? "--plp ${plp_prefix}" : "--sam $bam"

    """
    popscle demuxlet \\
        $input  \\
        --vcf ${donor_genotype} \\
        --out $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.best
    """
}
