process ANNOTSV_INSTALLANNOTATIONS {
    tag 'AnnotSV annotations'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2aa8892f5160501bd74eb44461ca6743dcf7132dd3b8bc769d09806585754b38/data' :
        'community.wave.seqera.io/library/annotsv:3.5.10--67b698fbe288f0ac' }"

    input:
    tuple val(meta), val(annotsv_version), val(exomiser_version)

    output:
    path "AnnotSV_annotations", emit: annotations
    tuple val("${task.process}"), val('annotsv'), eval("AnnotSV --version | sed 's/AnnotSV //'"), emit: versions_annotsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    annotsv_v = annotsv_version ?: ''
    exomiser_v = exomiser_version ?: ''

    if (( annotsv_v == '' && !(exomiser_v == '') ) || ( !(annotsv_v == '') && exomiser_v == '' ) ){
        error "Must supply both annotsv_version and exomiser_version or neither"
    }
    """
    INSTALL_annotations.sh  ${annotsv_v} ${exomiser_v}
    """

    stub:
    annotsv_v = annotsv_version ?: ''
    exomiser_v = exomiser_version ?: ''

    if (( annotsv_v == '' && !(exomiser_v == '') ) || ( !(annotsv_v == '') && exomiser_v == '' ) ){
        error "Must supply both annotsv_version and exomiser_version or neither"
    }
    """
    mkdir AnnotSV_annotations
    touch AnnotSV_annotations/stub_file.txt
    """
}
