process COOLER_BALANCE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.10.4--pyhdfd78af_0' :
        'biocontainers/cooler:0.10.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cool), val(resolution)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: cool
    tuple val("${task.process}"), val('cooler'), eval('cooler --version 2>&1 | sed "s/cooler, version //"'), emit: versions_cooler, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = resolution ? "::/resolutions/$resolution" : ""
    extension = cool.getExtension()
    if ("$cool" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    ln -s ${cool} ${prefix}.${extension}

    cooler balance \\
        $args \\
        -p ${task.cpus} \\
        ${prefix}.${extension}${suffix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = resolution ? "::/resolutions/$resolution" : ""
    extension = cool.getExtension()
    def creation_cmd = suffix.endsWith(".gz") ? "echo '' | gzip -c >" : "touch"
    """
    ${creation_cmd} ${prefix}.${extension}${suffix}
    """
}
