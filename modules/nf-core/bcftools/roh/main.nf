process BCFTOOLS_ROH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(af_file), path(af_file_tbi)
    path genetic_map
    path regions_file
    path samples_file
    path targets_file

    output:
    tuple val(meta), path("*.roh"), emit: roh
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def af_read   = af_file         ? "--AF-file ${af_file}"           : ''
    def gen_map   = genetic_map     ? "--genetic-map ${genetic_map}"   : ''
    def reg_file  = regions_file    ? "--regions-file ${regions_file}" : ''
    def samp_file = samples_file    ? "--samples-file ${samples_file}" : ''
    def targ_file = targets_file    ? "--targets-file ${targets_file}" : ''
    """
    bcftools \\
        roh \\
        $args \\
        $af_read \\
        $gen_map \\
        $reg_file \\
        $samp_file \\
        $targ_file \\
        -o ${prefix}.roh \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.roh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
