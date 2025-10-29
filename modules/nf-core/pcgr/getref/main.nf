process PCGR_GETREF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6d21ecd9fb81d5efb452bbcd06a72a23bb7dcc215610d1f92bf55dfe5a4eeee/data'
        : 'community.wave.seqera.io/library/coreutils_gzip_tar_wget:7fb7ade7a5b63d7a'}"

    input:
    tuple val(meta), val(bundleversion), val(genome)

    output:
    tuple val(meta), path("${bundleversion}"), emit: pcgrref
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bundle = "pcgr_ref_data.${bundleversion}.${genome}.tgz"
    """
    wget https://insilico.hpc.uio.no/pcgr/${bundle}
    gzip -dc ${bundle} | tar xvf -

    mkdir ${bundleversion}
    mv data/ ${bundleversion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    """
    mkdir ${bundleversion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}
