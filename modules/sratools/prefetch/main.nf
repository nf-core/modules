process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'
    label 'error_retry'

    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0' :
        'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0' }"

    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path("$id"), emit: sra
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = "/LIBS/GUID = \"${UUID.randomUUID().toString()}\"\\n/libs/cloud/report_instance_identity = \"true\"\\n"
    """
    eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"
    if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
        mkdir -p "\$(dirname "\${NCBI_SETTINGS}")"
        printf '${config}' > "\${NCBI_SETTINGS}"
    fi

    prefetch \\
        $args \\
        --progress \\
        $id

    vdb-validate $id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
