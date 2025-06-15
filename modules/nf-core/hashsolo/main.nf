process HASHSOLO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scanpy:1.11.2--e1f13faeaff34e5b':
        'community.wave.seqera.io/library/scanpy:1.11.2--cc70dffaa98ddb98' }"

    input:
    tuple val(meta), path(input_h5ad)
    val params

    output:
    tuple val(meta), path("test/*.h5ad"), emit: h5ad
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('hashsolo.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p test
    touch test/${prefix}_hashsolo.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: "STUB"
    END_VERSIONS
    """
}
