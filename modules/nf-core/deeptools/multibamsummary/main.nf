process DEEPTOOLS_MULTIBAMSUMMARY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bams), path(bais), val(labels)

    output:
    tuple val(meta), path("*.npz") , emit: matrix
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "all_bam"
    def label  = labels ? "--labels ${labels.join(' ')}" : ''
    """
    multiBamSummary bins \\
        $args \\
        $label \\
        --bamfiles ${bams.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName ${prefix}.bamSummary.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBamSummary --version | sed -e "s/multiBamSummary //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_bam"
    """
    touch ${prefix}.bamSummary.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBamSummary --version | sed -e "s/multiBamSummary //g")
    END_VERSIONS
    """
}
