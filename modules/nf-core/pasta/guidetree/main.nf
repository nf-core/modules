process PASTA_GUIDETREE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-968f27e3bbe4396c3ca86df137949c3cd3b69a97:83d0b2752dbcce91d17a4528e065feb7ae01ee25':
        'biocontainers/mulled-v2-968f27e3bbe4396c3ca86df137949c3cd3b69a97:83d0b2752dbcce91d17a4528e065feb7ae01ee25' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tre"), emit: tree
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "output": "${meta.id}"
    """
    run_pasta.py \\
        --num-cpus $task.cpus \\
        -i <(pigz -cdf $fasta) \\
        -o ./ \\
        -j $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PASTA: \$(run_pasta.py --version | sed 's/^PASTA\\wv//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" > ${prefix}.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PASTA: \$(run_pasta.py --version | sed 's/^PASTA\\wv//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
