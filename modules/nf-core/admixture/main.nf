// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join


// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process ADMIXTURE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::admixture=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/admixture:1.3.0--0':
        'quay.io/biocontainers/admixture:1.3.0--0' }"

    input:
    tuple val(meta), path(input_file)
    val(K)

    output:
    tuple val(meta), path("*.Q"), emit: Q-ancestry-fractions
    tuple val(meta), path("*.P"), emit: P-allele-frequencies
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    admixture \\
        $input_file \\
        $bed \\
        $K \\
        -J $task.cpus \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(echo \$(admixture 2>&1) | head -n 1 |  sed -n 's/.*Version \([0-9.]*\).*/\1/p' )
    END_VERSIONS
    """
}
