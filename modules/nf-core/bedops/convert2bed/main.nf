process BEDOPS_CONVERT2BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h4ac6f70_2':
        'biocontainers/bedops:2.4.41--h4ac6f70_2' }"

    input:
    tuple val(meta), path(in_file)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = in_file.getExtension()
    """
    convert2bed \\
        $args \\
        -i $format \\
        < $in_file \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops/convert2bed: \$(convert2bed --version | grep vers | sed 's/^.*.version:  //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops/convert2bed: \$(convert2bed --version | grep vers | sed 's/^.*.version:  //')
    END_VERSIONS
    """
}
