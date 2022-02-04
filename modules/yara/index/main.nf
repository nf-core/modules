process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::yara=1.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yara:1.0.2--2' :
        'quay.io/biocontainers/yara:1.0.2--2' }"

    input:
    path fasta

    output:
    path "yara"        , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir yara

    yara_indexer \\
        $fasta \\
        -o "yara"

    mv *.{lf,rid,sa,txt}.* yara
    cp $fasta yara/yara.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yara: \$(echo \$(yara_indexer --version 2>&1) | sed 's/^.*yara_indexer version: //; s/ .*\$//')
    END_VERSIONS
    """
}
