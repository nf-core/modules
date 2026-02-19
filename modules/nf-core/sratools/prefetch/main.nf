process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.1.0--h9f5acd7_0' :
        'biocontainers/sra-tools:3.1.0--h9f5acd7_0' }"

    input:
    tuple val(meta), val(id)
    path ncbi_settings
    path certificate

    output:
    tuple val(meta), path("${id}", type: 'dir'), emit: sra
    tuple val("${task.process}"), val('sratools'), eval("prefetch --version 2>&1 | grep -Eo '[0-9.]+'"), topic: versions, emit: versions_sratools
    tuple val("${task.process}"), val('curl'), eval("curl --version | head -n 1 | sed 's/^curl //; s/ .*\$//'"), topic: versions, emit: versions_curl

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    args2 = task.ext.args2 ?: '5 1 100'  // <num retries> <base delay in seconds> <max delay in seconds>
    def cert_arg = certificate ?
        (certificate.name.endsWith('.jwt') ? "--perm ${certificate}" :
        certificate.name.endsWith('.ngc') ? "--ngc ${certificate}" : '') : ''
    def final_args = "${args} ${cert_arg}".trim()

    def retry_script_content = file("${moduleDir}/templates/retry_with_backoff.sh").text

    """
    export NCBI_SETTINGS="\${PWD}/${ncbi_settings}"
    export SRA_ID="${id}"
    export PROCESS_NAME="${task.process}"
    export PREFETCH_ARGS="${final_args}"
    export RETRY_ARGS="${args2}"

    cat > retry_with_backoff.sh <<'SCRIPT_EOF'
${retry_script_content}
SCRIPT_EOF

    bash retry_with_backoff.sh
    """

    stub:
    """
    mkdir ${id}
    touch ${id}/${id}.sra
    """
}
