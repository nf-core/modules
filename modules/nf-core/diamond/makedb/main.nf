process DIAMOND_MAKEDB {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::diamond=2.0.15" : null)
    def container_image = "diamond:2.0.15--hb97b32f_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    path fasta

    output:
    path "${fasta}.dmnd", emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    diamond \\
        makedb \\
        --threads $task.cpus \\
        --in  $fasta \\
        -d $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
