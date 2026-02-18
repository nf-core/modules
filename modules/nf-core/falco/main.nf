process FALCO {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/falco:1.2.5--h077b44d_0'
        : 'biocontainers/falco:1.2.5--h077b44d_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('falco'), eval("falco --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_falco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (reads.toList().size() == 1) {
        """
        falco ${args} --threads ${task.cpus} ${reads} -D ${prefix}_fastqc_data.txt -S ${prefix}_summary.txt -R ${prefix}_report.html
        """
    }
    else {
        """
        falco ${args} --threads ${task.cpus} ${reads}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_data.txt
    touch ${prefix}_fastqc_data.html
    touch ${prefix}_summary.txt
    """
}
