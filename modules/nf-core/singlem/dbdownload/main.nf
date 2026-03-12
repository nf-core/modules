process SINGLEM_DBDOWNLOAD {
    tag '$bam'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    
    path dna_sequence

    output:
    path taxonomy, emit: taxonomy
    tuple val("${task.process}"), val('singlem'), eval("singlem --version"), topic: versions, emit: versions_singlem

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    singlem \\
        $args \\
        -@ $task.cpus \\
        $dna_sequence \\
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    echo $args
    
    touch ${prefix}.fastq
    touch ${prefix}.fasta
    """
}
