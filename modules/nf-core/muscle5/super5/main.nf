process MUSCLE5_SUPER5 {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8eb01a3c2755c935d070dd03ff2dee698eeb4466:ceb6e65e00346ed20d0d8078dddf9858a7af0fe2-0':
        'biocontainers/mulled-v2-8eb01a3c2755c935d070dd03ff2dee698eeb4466:ceb6e65e00346ed20d0d8078dddf9858a7af0fe2-0' }"

    input:
    tuple val(meta), path(fasta)
    val(compress)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    prefix = args.contains('-perm all') ? "${prefix}@" : "${prefix}"
    """
    muscle \\
        -super5 ${fasta} \\
        -output ${prefix}.aln \\
        ${args} \\
        -threads ${task.cpus}

    # output may be multiple files if -perm all is set
    for f in *.aln; do
        pigz -p ${task.cpus} \$f
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muscle: \$(muscle -version | head -n 1 | cut -d ' ' -f 2 | sed 's/.linux64//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muscle: \$(muscle -version | head -n 1 | cut -d ' ' -f 2 | sed 's/.linux64//')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
