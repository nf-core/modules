// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.


process DEEPTOOLS_MULTIBAMSUMMARY {
    tag "$meta.id"
    label 'process_high' 

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bams), path(bais), val(labels)
    
    output:
    tuple val(meta), path("*.npz") , emit: matrix
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiBamSummary bins \\
        $args \\
        --bamfiles ${bams.join(' ')} \\
        --labels ${labels.join(' ')} \\
        --numberOfProcessors $task.cpus \\
        --outFileName all_bam.bamSummary.npz \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBamSummary --version | sed -e "s/multiBamSummary //g")
    END_VERSIONS
    """
}