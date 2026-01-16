process SAMTOOLS_QUICKCHECK {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'terminate'   // Nextflow default; but explicit here

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0':
        'biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    /*  There is not generated output files.
        The sole purpose of 'samtools quickcheck' is to return an error if the mapping file is malformatted.
        The purpose in of the Nextflow process is to break the main workflow (by default) if the script returns a non-zero exit code.
        The nf-core components guidelines (https://nf-co.re/docs/guidelines/components/modules) do not specify that there must be output files.
    */
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools quickcheck \\
        $args \\
        $bam
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo $args
    """
}
