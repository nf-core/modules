process RASUSA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::rasusa=0.3.0" : null)
    def container_image = "/rasusa:0.3.0--h779adbc_1"
                                         container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(reads), val(genome_size)
    val   depth_cutoff

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output   = meta.single_end ? "--output ${prefix}.fastq.gz" : "--output ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz"
    """
    rasusa \\
        $args \\
        --coverage $depth_cutoff \\
        --genome-size $genome_size \\
        --input $reads \\
        $output
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version 2>&1 | sed -e "s/rasusa //g")
    END_VERSIONS
    """
}
