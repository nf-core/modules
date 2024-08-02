process PASTA_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-968f27e3bbe4396c3ca86df137949c3cd3b69a97:83d0b2752dbcce91d17a4528e065feb7ae01ee25':
        'biocontainers/mulled-v2-968f27e3bbe4396c3ca86df137949c3cd3b69a97:83d0b2752dbcce91d17a4528e065feb7ae01ee25' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(tree)
    val(compress)

    output:
    tuple val(meta), path("$prefix.aln{.gz,}"), emit: alignment
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def loadtree = tree ? "-t $tree" : ''
    def compress_output = compress ? "pigz -p ${task.cpus} ${prefix}.aln" : ""
    // pasta opens and writes multiple output files multiple times; compress after the fact
    """
    run_pasta.py \\
        --num-cpus $task.cpus \\
        -i <(pigz -cd $fasta) \\
        -o ./ \\
        -j $prefix \\
        --alignment-suffix aln\\
        $loadtree \\
        $args

    # pasta writes many temp files into the dir, compress only output one
    $compress_output

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
    echo "" | gzip > ${prefix}.aln${compress ? '.gz' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PASTA: \$(run_pasta.py --version | sed 's/^PASTA\\wv//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
