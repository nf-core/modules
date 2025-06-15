process HASHSOLO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0':
        'biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_h5ad)
    val params

    output:
    tuple val(meta), path("test/*.h5ad"), emit: h5ad
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    template('hashsolo.py')

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p test
    touch test/${prefix}_hashsolo.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: "STUB"
    END_VERSIONS
    """
}
