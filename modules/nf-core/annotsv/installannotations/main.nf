process ANNOTSV_INSTALLANNOTATIONS {
    tag 'AnnotSV annotations'
    label 'process_single'

    conda "bioconda::annotsv=3.3.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annotsv:3.3.6--py311hdfd78af_0' :
        'biocontainers/annotsv:3.3.6--py311hdfd78af_0' }"

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV --version | sed 's/AnnotSV //'))
    END_VERSIONS
    """
}
