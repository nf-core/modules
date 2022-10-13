process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :

    input:
    path fasta

    output:
    path "bwa"         , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${fasta.baseName} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
