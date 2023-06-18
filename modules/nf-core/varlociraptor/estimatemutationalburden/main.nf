process VARLOCIRAPTOR_ESTIMATEMUTATIONALBURDEN {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::varlociraptor=8.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor:8.1.1--hc349b7f_0':
        'biocontainers/varlociraptor:8.1.1--hc349b7f_0' }"

    input:
    tuple val(meta), path(vcf)
    val(plot_mode)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    varlociraptor estimate mutational-burden \\
        --plot-mode ${plot_mode} \\
        $args \\
        < ${vcf} \\
        > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
