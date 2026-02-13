process DIANN_INSILICOLIBRARYGENERATION {
    tag "$fasta.Name"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/diann/v1.8.1_cv1/diann_v1.8.1_cv1.img' :
        'docker.io/biocontainers/diann:v1.8.1_cv1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.predicted.speclib"), emit: predict_speclib
    tuple val(meta), path("*.log"),               emit: log
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN_INSILICOLIBRARYGENERATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    diann \\
            --fasta ${fasta} \\
            --fasta-search \\
            --threads ${task.cpus} \\
            --predictor \\
            --gen-spec-lib \\
            ${args}

    cp *lib.log.txt ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diann: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN_INSILICOLIBRARYGENERATION module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.predicted.speclib
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diann: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """
}
