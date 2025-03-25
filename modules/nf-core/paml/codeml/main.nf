process PAML_CODEML {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paml:4.10.7--h7b50bb2_2':
        'quay.io/biocontainers/paml:4.10.7--h7b50bb2_2' }"

    input:
    tuple val(meta), path(phy)
    path tree
    path ctl

    output:
    tuple val(meta), path("result.txt")  , emit: out_txt
    path "out.log"                       , emit: log
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat  ${ctl} | awk '/outfile/ { print "outfile = result.txt"; next }{ print }' | tee ${ctl} 2>&1

    codeml "${ctl}" | tee out.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paml: \$(head out.log | grep 'CODONML')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch result.txt
    echo "CODONML in paml version 4.10.7, June 2023" > out.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paml: \$(head out.log | grep 'CODONML')
    END_VERSIONS
    """
}
