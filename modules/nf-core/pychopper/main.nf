process PYCHOPPER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pychopper:2.7.10--pyhdfd78af_0':
        'quay.io/biocontainers/pychopper:2.7.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.out.fastq.gz"), emit: fastq
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    tuple val("${task.process}"), val('pychopper'), val('2.7.10'), emit: versions_pychopper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pychopper \\
        ${args} \\
        -t ${task.cpus} \\
        ${fastq} \\
        ${prefix}.out.fastq

    gzip -f ${prefix}.out.fastq > ${prefix}.out.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.out.fastq.gz
    """
}
