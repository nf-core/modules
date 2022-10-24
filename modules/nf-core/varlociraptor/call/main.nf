process VARLOCIRAPTOR_CALL {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::varlociraptor=5.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor5.3.3--hc349b7f_0':
        'quay.io/biocontainers/varlociraptor5.3.3--hc349b7f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //'  ))
    END_VERSIONS
    """
}
