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
    tuple val(meta), path(id, type: 'dir'), emit: sra
    path 'versions.yml'      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args = task.ext.args ?: ''
    args2 = task.ext.args2 ?: '5 1 100'  // <num retries> <base delay in seconds> <max delay in seconds>
    if (certificate) {
        if (certificate.toString().endsWith('.jwt')) {
            args += " --perm ${certificate}"
        }
        else if (certificate.toString().endsWith('.ngc')) {
            args += " --ngc ${certificate}"
        }
    }

    template 'retry_with_backoff.sh'

    stub:
    """
    mkdir $id
    touch $id/${id}.sra

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
        curl: \$(curl --version | head -n 1 | sed 's/^curl //; s/ .*\$//')
    END_VERSIONS
    """
}
