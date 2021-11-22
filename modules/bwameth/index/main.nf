process BWAMETH_INDEX {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwameth=0.2.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bwameth:0.2.2--py_1"
    } else {
        container "quay.io/biocontainers/bwameth:0.2.2--py_1"
    }

    input:
    path fasta, stageAs: "bwameth/*"

    output:
    path "bwameth"      , emit: index
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bwameth.py index $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bwameth.py --version 2>&1) | cut -f2 -d" ")
    END_VERSIONS
    """
}
