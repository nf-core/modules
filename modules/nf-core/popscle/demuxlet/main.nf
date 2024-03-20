process POPSCLE_DEMUXLET {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/popscle:0.1beta--h2c78cec_0' :
        'biocontainers/popscle:0.1beta--h2c78cec_0' }"

    input:
    tuple val(meta), path(plp_input), path(bam), path(donor_genotype)

    output:
    tuple val(meta), path('*.best'), emit: demuxlet_result
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = plp_input ? "--plp ${plp_input}" : "--sam $bam"
    def VERSION ='0.1'
    """
    popscle demuxlet \\
        $input  \\
        --vcf ${donor_genotype} \\
        --out $prefix \\
	    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle demuxlet: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = plp_input ? "--plp ${plp_input}" : "--sam $bam"
    def VERSION ='0.1'
    """
    touch ${prefix}.best

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        popscle demuxlet: $VERSION
    END_VERSIONS
    """
}
