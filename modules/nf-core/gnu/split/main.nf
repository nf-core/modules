process GNU_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path( "*split*" ), emit: split
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = input.extension
    if (suffix == 'gz') {
        def next_suffix = file(input.baseName).getExtension()
        """
        gunzip -c ${input} | split ${args} --additional-suffix=.${next_suffix} - ${prefix}.split.
        gzip ${prefix}.split.*

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gnu: \$(split --version |& sed '1!d ; s/split (GNU coreutils) //')
        END_VERSIONS
        """
    } else {
        """
        split ${args} --additional-suffix=.${suffix} ${input} ${prefix}.split.

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gnu: \$(split --version |& sed '1!d ; s/split (GNU coreutils) //')
        END_VERSIONS
        """
    }

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.split.000.csv ${prefix}.split.001.csv ${prefix}.split.002.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gnu: \$(split --version |& sed '1!d ; s/split (GNU coreutils) //')
    END_VERSIONS
    """
}
