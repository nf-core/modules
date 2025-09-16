process ANNOTSV_INSTALLANNOTATIONS {
    tag 'AnnotSV annotations'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/363f212881f1b2f5c3395a6c7d1270694392e3a6f886e46e091e83527fed9b6b/data' :
        'community.wave.seqera.io/library/annotsv:3.5.3--71a461cb86d570b7' }"

    output:
    path "AnnotSV_annotations", emit: annotations
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    TAG=\$(AnnotSV --version 2>/dev/null | sed 's/AnnotSV //' || echo "master")
    INSTALL_annotations.sh \$TAG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV --version | sed 's/AnnotSV //'))
    END_VERSIONS
    """

    stub:
    """
    mkdir AnnotSV_annotations
    touch AnnotSV_annotations/stub_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV --version | sed 's/AnnotSV //'))
    END_VERSIONS
    """
}
