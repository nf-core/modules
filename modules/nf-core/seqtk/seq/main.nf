process SEQTK_SEQ {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
        def container_image = "/seqtk:1.3--h5bf99c6_3"
                                            container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.gz")     , emit: fastx
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz/ || "$args" ==~ /\-[aA]/ ) {
        extension = "fasta"
    }
    """
    seqtk \\
        seq \\
        $args \\
        $fastx | \\
        gzip -c > ${prefix}.seqtk-seq.${extension}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
