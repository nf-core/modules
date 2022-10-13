process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwameth=0.2.2" : null)
    def container_image = "/bwameth:0.2.2--py_1"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
