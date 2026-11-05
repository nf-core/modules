// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ANNOTSV_INSTALLANNOTATIONS {
    tag 'AnnotSV annotations'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/664c800960f1b08bcae320959b20db6b606f3c5de2488d63376eed58f32b1a10/data'
        : 'community.wave.seqera.io/library/annotsv:3.5.10--28ac84a07bcbf122'}"

    input:
    tuple val(meta), val(annotsv_version), val(exomiser_version)

    output:
    tuple val(meta), path("AnnotSV_annotations"), emit: annotations
    tuple val("${task.process}"), val('annotsv'), eval("AnnotSV --version | sed 's/AnnotSV //'"), emit: versions_annotsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    annotsv_v = annotsv_version ?: ''
    exomiser_v = exomiser_version ?: ''

    if ( (!annotsv_v) != (!exomiser_v) ) {
      error "Must supply both annotsv_version and exomiser_version or neither"
  }
    """
    INSTALL_annotations.sh  ${annotsv_v} ${exomiser_v}
    """

    stub:
    annotsv_v = annotsv_version ?: ''
    exomiser_v = exomiser_version ?: ''

    if ( (!annotsv_v) != (!exomiser_v) ) {
      error "Must supply both annotsv_version and exomiser_version or neither"
  }
    """
    mkdir AnnotSV_annotations
    touch AnnotSV_annotations/stub_file.txt
    """
}
