process NANOLYSE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::nanolyse=1.2.0" : null)
    def container_image = "nanolyse:1.2.0--py_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


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
