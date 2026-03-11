process BCFTOOLS_ROH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(af_file), path(af_file_tbi)
    path genetic_map
    path regions_file
    path samples_file
    path targets_file

    output:
    tuple val(meta), path("*.roh"), emit: roh
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def af_read = af_file ? "--AF-file ${af_file}" : ''
    def gen_map = genetic_map ? "--genetic-map ${genetic_map}" : ''
    def reg_file = regions_file ? "--regions-file ${regions_file}" : ''
    def samp_file = samples_file ? "--samples-file ${samples_file}" : ''
    def targ_file = targets_file ? "--targets-file ${targets_file}" : ''
    """
    bcftools \\
        roh \\
        ${args} \\
        ${af_read} \\
        ${gen_map} \\
        ${reg_file} \\
        ${samp_file} \\
        ${targ_file} \\
        -o ${prefix}.roh \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.roh
    """
}
