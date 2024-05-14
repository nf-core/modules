// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process RTN_TNI {
    debug true
    tag "{$expression_matrix.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtn:2.26.0--r43hdfd78af_0':
        'biocontainers/bioconductor-rtn:2.26.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(expression_matrix)

    output:
    path "tni.rds"                     , emit: tni
    path "tni_permutated.rds"          , emit: tni_perm
    path "tni_bootstrapped.rds"        , emit: tni_bootstrap
    path "tni_filtered.rds"            , emit: tni_filtered
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    template 'rtn_tni.r'

    stub:
    def args = task.ext.args ?: ''

    """
    touch tni.rds
    touch tni_permutated.rds
    touch tni_bootstrapped.rds
    touch tni_filtered.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-rtn: \$(Rscript -e "suppressWarnings(library(RTN)); cat(as.character(packageVersion('RTN')))")
    END_VERSIONS
    """
}
