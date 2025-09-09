process XZ_DECOMPRESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0' :
        'biocontainers/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("$decompressed_file"), emit: file
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    decompressed_file = archive.toString().replaceAll(".xz\$", "")
    """
    # Note 1: needs --stdout for xz --decompress to avoid two issues:
    #   1. xz: ${archive}: Is a symbolic link, skipping
    #   2. xz: ${archive}: Cannot set the file group: Operation not permitted
    # Note 2: using several threads in xz --decompress will only work on files that contain multiple blocks with size
    # information in block headers.  All files compressed in multi-threaded mode meet this condition.
    xz -T ${task.cpus} --decompress --stdout ${args} ${archive} > ${decompressed_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xz: \$(xz --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    decompressed_file = archive.toString().replaceAll(".xz\$", "")
    """
    touch "${decompressed_file}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xz: \$(xz --version | head -n1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
