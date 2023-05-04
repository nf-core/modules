process NANOLYSE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::nanolyse=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanolyse:1.2.0--py_0' :
        'biocontainers/nanolyse:1.2.0--py_0' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "*.log"                       , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c $fastq | NanoLyse -r $fasta | gzip > ${prefix}.fastq.gz
    mv NanoLyse.log ${prefix}.nanolyse.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanolyse: \$(NanoLyse --version 2>&1 | sed -e "s/NanoLyse //g")
    END_VERSIONS
    """
}
