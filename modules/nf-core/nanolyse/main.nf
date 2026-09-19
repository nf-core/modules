process NANOLYSE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanolyse:1.2.0--py_0' :
        'quay.io/biocontainers/nanolyse:1.2.0--py_0' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.log"                       , emit: log
    tuple val("${task.process}"), val('nanolyse'), eval('NanoLyse --version 2>&1 | sed -e "s/NanoLyse //g"'), emit: versions_nanolyse, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c $fastq | NanoLyse -r $fasta | gzip > ${prefix}.fastq.gz
    mv NanoLyse.log ${prefix}.nanolyse.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    touch ${prefix}.nanolyse.log
    """
}
