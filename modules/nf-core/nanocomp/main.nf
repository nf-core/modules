process NANOCOMP {
    label 'process_medium'

    conda "bioconda:nanocomp=1.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp:1.21.0--pyhdfd78af_0':
        'quay.io/biocontainers/nanocomp:1.21.0--pyhdfd78af_0' }"

    input:
    val filetype
    path filelist

    output:
    path "*.html", emit: html_nanocomp_output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if ! [[ $filetype = @(fastq|fasta|summary|bam) ]]; then echo "Input should only be fastq, fasta, summary or bam type, but is $filetype" && exit 1; fi
    NanoComp \\
        --$filetype $filelist \\
        --verbose \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(echo \$(NanoComp --version 2>&1) | sed 's/^.*NanoComp //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
