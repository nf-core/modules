process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::bwameth=0.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwameth:0.2.2--py_1' :
        'biocontainers/bwameth:0.2.2--py_1' }"

    input:
    path fasta, stageAs: "bwameth/*"

    output:
    path "bwameth"      , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bwameth.py index $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(echo \$(bwameth.py --version 2>&1) | cut -f2 -d" ")
    END_VERSIONS
    """
}
