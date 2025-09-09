process PYDAMAGE_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pydamage:1.0--pyhdfd78af_0' :
        'biocontainers/pydamage:1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}_pydamage_results.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=./tmp
    export MPLCONFIGDIR=./tmp

    pydamage \\
        analyze \\
        $args \\
        -p $task.cpus \\
        $bam

    mv pydamage_results/pydamage_results.csv ${prefix}_pydamage_results.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pydamage: \$(pydamage --version | sed -n 's/pydamage, version \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=./tmp # exports required for versions to resolve correctly
    export MPLCONFIGDIR=./tmp

    touch ${prefix}_pydamage_results.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pydamage: \$(echo \$(pydamage --version 2>&1) | sed -e 's/pydamage, version //g')
    END_VERSIONS
    """

}
