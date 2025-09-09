process ANNOTSV_INSTALLANNOTATIONS {
    tag 'AnnotSV annotations'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/df/df5e87ce610d595afe5f1e4017c255e53590cfa9b156c43310d637b36bfd12b3/data' :
        'community.wave.seqera.io/library/annotsv:3.4.6--c5a6da1bf4c2c8e5' }"

    output:
    path "AnnotSV_annotations", emit: annotations
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    INSTALL_annotations.sh

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
