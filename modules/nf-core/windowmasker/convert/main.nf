process WINDOWMASKER_CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0':
        'biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${output}"), emit: converted
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfmt =    args.contains('-sformat binary')  ? 'binary'  :
                    args.contains('-sformat oascii')  ? 'oascii'  :
                    args.contains('-sformat obinary') ? 'obinary' :
                    'ascii'
    output  = "${prefix}.${outfmt}"
    """
    windowmasker -convert \\
        -in $counts \\
        -out $output \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfmt =    args.contains('-sformat binary')  ? 'binary'  :
                    args.contains('-sformat oascii')  ? 'oascii'  :
                    args.contains('-sformat obinary') ? 'obinary' :
                    'ascii'
    output  = "${prefix}.${outfmt}"
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """
}
