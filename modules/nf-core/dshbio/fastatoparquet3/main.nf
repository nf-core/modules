process DSHBIO_FASTATOPARQUET3 {
    tag '$fasta'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dsh-bio:2.4--hdfd78af_0':
        'biocontainers/dsh-bio:2.4--hdfd78af_0' }"

    input:
    path fasta

    output:
    path "*.sequences.parquet", emit: parquet
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.getBaseName(2)}"

    """
    dsh-bio \\
        fasta-to-parquet3 \\
        $args \\
        -i $fasta \\
        -o ${prefix}.sequences.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dshbio: \$(dsh-bio --version 2>&1 | grep -o 'dsh-bio-tools .*' | cut -f2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.getBaseName(2)}"

    """
    mkdir -p ${prefix}.sequences.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dshbio: \$(dsh-bio --version 2>&1 | grep -o 'dsh-bio-tools .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
