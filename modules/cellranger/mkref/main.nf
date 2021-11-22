process CELLRANGER_MKREF {
    tag 'mkref'
    label 'process_high'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "nfcore/cellranger:6.0.2"

    input:
    path fasta
    path gtf
    val(reference_name)

    output:
    path "versions.yml"     , emit: versions
    path "${reference_name}", emit: reference

    script:
    def args = task.ext.args ?: ''
    """
    cellranger mkref \\
        --genome=${reference_name} \\
        --fasta=${fasta} \\
        --genes=${gtf}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    }
    container "nfcore/cellranger:6.0.2"

    input:
    path fasta
    path gtf
    val(reference_name)

    output:
    path "versions.yml"     , emit: versions
    path "${reference_name}", emit: reference

    script:
    def args = task.ext.args ?: ''
    """
    cellranger mkref \\
        --genome=${reference_name} \\
        --fasta=${fasta} \\
        --genes=${gtf}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '{' :
        'process_high' }"
}
