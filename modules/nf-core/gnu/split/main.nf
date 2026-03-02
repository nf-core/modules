process GNU_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5':
        'biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path( "${prefix}.*" )  , emit: split
    tuple val("${task.process}"), val('coreutils'), eval("sort --version |& sed '1!d ; s/sort (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}.split"
    def suffix  = input.extension
    if (suffix == 'gz') {
        def next_suffix = file(input.baseName).getExtension()
        """
        gunzip -c ${input} | split ${args} --additional-suffix=.${next_suffix} - ${prefix}.
        gzip ${prefix}.*

        """
    } else {
        """
        split ${args} --additional-suffix=.${suffix} ${input} ${prefix}.

        """
    }

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}.split"
    """
    touch ${prefix}.000.csv ${prefix}.001.csv ${prefix}.002.csv

    """
}
