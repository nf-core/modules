process RSEQC_BAMSTAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'" : null)
    def container_image = "/rseqc:3.0.1--py37h516909a_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: txt
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_stat.py \\
        -i $bam \\
        $args \\
        > ${prefix}.bam_stat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(bam_stat.py --version | sed -e "s/bam_stat.py //g")
    END_VERSIONS
    """
}
