process CUSTOM_CATADDITIONALFASTA {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta), path(gtf)
    tuple val(meta2), path(add_fasta)
    val(biotype)

    output:
    tuple val(meta), path("*/*.fasta"), emit: fasta
    tuple val(meta), path("*/*.gtf")  , emit: gtf
    path "versions.yml"               , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix2 = task.ext.prefix2 ?: "${meta2.id}"
    prefix = task.ext.prefix ?: "${meta.id}_${prefix2}"

    """
    echo $prefix2
    echo $prefix
    """

    template 'fasta2gtf.py'

    stub:
    prefix2 = task.ext.prefix2 ?: "${meta2.id}"
    prefix = task.ext.prefix ?: "${meta.id}_${prefix2}"
    """
    mkdir out
    touch out/${prefix}.fasta
    touch out/${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
