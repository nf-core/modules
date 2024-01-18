process LAI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ltr_retriever:2.9.0--hdfd78af_2':
        'biocontainers/ltr_retriever:2.9.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    path pass_list
    path annotation_out
    path monoploid_seqs

    output:
    tuple val(meta), path("*.LAI.log")  , emit: log
    tuple val(meta), path("*.LAI.out")  , emit: lai_out     , optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args     ?: ''
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    def monoploid_param = monoploid_seqs    ? "-mono $monoploid_seqs"                       : ''
    def lai_output_name = monoploid_seqs    ? "${annotation_out}.${monoploid_seqs}.out.LAI" : "${annotation_out}.LAI"
    def VERSION         = 'beta3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # Remove comments from genome fasta,
    # otherwise LAI triggers its sequence name change logic

    sed \\
        '/^>/ s/\\s.*\$//' \\
        $fasta \\
        > for_lai_no_comments.fsa

    LAI \\
        -genome for_lai_no_comments.fsa \\
        -intact $pass_list \\
        -all $annotation_out \\
        -t $task.cpus \\
        $monoploid_param \\
        $args \\
        > "${prefix}.LAI.log"

    mv \\
        $lai_output_name \\
        "${prefix}.LAI.out" \\
        || echo "LAI did not produce the output file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lai: $VERSION
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = 'beta3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.LAI.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lai: $VERSION
    END_VERSIONS
    """
}
