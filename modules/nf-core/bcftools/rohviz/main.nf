process BCFTOOLS_ROHVIZ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/2397e92afdb695a248544593b5dcd19d8a0321dd3fd7f85cacd80de6cd56290d/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_less:6e82455ad302b6f7' }"

    input:
    tuple val(meta), path(roh)
    tuple val(meta1), path(vcf)
    path regions_list
    path samples_file

    output:
    tuple val(meta), path("*.html") , emit: html
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools
    tuple val("${task.process}"), val('less'), eval("less --version | head -n1 | sed 's/less //'"), topic: versions, emit: versions_less

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def regions   = regions_list    ? "--regions ${regions_list}"      : ''
    def samp_file = samples_file    ? "--samples-file ${samples_file}" : ''

    // Set TERM=dumb to prevent zless/less from failing in non-interactive container environments

    """
    env TERM=dumb roh-viz \\
        $args \\
        -i ${roh} \\
        -v ${vcf}\\
        $regions \\
        $samp_file \\
        -o ${prefix}.roh-viz.html
    """

    stub:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "args: ${args}"
    touch ${prefix}.roh-viz.html
    """
}
