process GGET_GGET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7ca1a63d986bc5b68d0543ce326a703265d4a3645be7e05f9af94f2992aa8352/data':
        'community.wave.seqera.io/library/gget:0.29.1--3b5a50589bc0feb3' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*[!versions.yml][!${prefix}.${extension}]*"), emit: files , optional: true
    tuple val(meta), path("${prefix}.${extension}")                    , emit: output, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inputs = files ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.startsWith("pdb")  ? "pdb" :
                args.contains("-csv")   ? "csv" :
                "json"
    """
    export MPLCONFIGDIR=\$PWD/.tmp
    gget \\
        $args \\
        -o ${prefix}.${extension} \\
        $inputs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: \$(echo \$(gget --version 2>&1 | sed 's/gget version: //g'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.startsWith("pdb")  ? "pdb" :
                args.contains("-csv")   ? "csv" :
                "json"
    """
    export MPLCONFIGDIR=\$PWD/.tmp    #included in stub to stop errors in the version string
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: \$(echo \$(gget --version 2>&1 | sed 's/gget version: //g'))
    END_VERSIONS
    """
}
