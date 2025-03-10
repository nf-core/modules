process SAMTOOLS_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta),  path(input_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional:true, emit: cram
    path  "versions.yml"                                  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    samtools \\
        cat \\
        $args \\
        -o ${prefix}.${file_type} \\
        $input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    """
    touch ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
