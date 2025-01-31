process FIND_UNPIGZ {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fd226561e12b32bcacdf4f5ff74577e76233adf52ae5cbc499a2cdfe0e27d82/data'
        : 'community.wave.seqera.io/library/findutils_pigz:c4dd5edc44402661'}"

    input:
    tuple val(meta), path(files_in, stageAs: 'gzipped/*', arity: '1..*')

    output:
    tuple val(meta), path("${prefix}.*"), emit: file_out
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    if (files_in.any { file -> !file.name.endsWith('.gz') }) {
        error("All files provided to this module must be gzipped (and have the .gz extension).")
    }

    """
    while IFS= read -r -d \$'\\0' file; do
        unpigz \\
            ${args} \\
            -cd \\
            --processes ${task.cpus} \\
            \$file \\
            > ${prefix}.\$( basename \$file .gz )
    done < <( find gzipped/ -name '*.gz' -print0 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find: \$( find --version | sed '1!d; s/.* //' )
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.test_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find: \$( find --version | sed '1!d; s/.* //' )
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
