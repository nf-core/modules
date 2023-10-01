process EIDO_VALIDATE {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::eido=0.1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/eido/0.1.9_cv2/eido_0.1.9_cv2.sif' :
        'docker.io/biocontainers/eido:0.1.9_cv2' }"

    input:
    path samplesheet
    path schema
    path pep_input_base_dir

    output:
    path "versions.yml"  , emit: versions
    path "*.log"         , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "validation"
    """
    eido validate $args $samplesheet -s $schema -e > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eido: \$(echo \$(eido --version 2>&1) | sed 's/^.*eido //;s/ .*//' )
    END_VERSIONS
    """
}
