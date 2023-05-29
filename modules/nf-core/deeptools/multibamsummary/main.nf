process DEEPTOOLS_MULTIBAMSUMMARY {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bams), path(bais), val(labels)

    output:
    tuple val(meta), path("*.npz") , emit: matrix
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def label = labels ? "--labels ${labels.join(' ')}" : ''
    """
    multiBamSummary bins \\
        $args \\
        $label \\
        --bamfiles ${bams.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName all_bam.bamSummary.npz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBamSummary --version | sed -e "s/multiBamSummary //g")
    END_VERSIONS
    """
}
