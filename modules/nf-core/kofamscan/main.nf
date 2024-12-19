process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
        'biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    path profiles
    path ko_list

    output:
    tuple val(meta), path("*.txt"), optional: true, emit: txt
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format detail-tsv") ? "tsv" :
                    "txt"
    """
    exec_annotation \\
        -p $profiles \\
        -k $ko_list \\
        $args \\
        --cpu $task.cpus \\
        -o ${prefix}.${extension} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kofamscan: \$(echo \$(exec_annotation --version 2>&1) | sed 's/^.*exec_annotation //;')
    END_VERSIONS
    """
}
