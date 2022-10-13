process HMMER_HMMALIGN {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    def container_image = "/hmmer:3.3.2--h1b792b2_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(fasta)
    path hmm

    output:
    tuple val(meta), path("*.sthlm.gz"), emit: sthlm
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hmmalign \\
        $args \\
        $hmm \\
        $fasta | gzip -c > ${prefix}.sthlm.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmalign -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
