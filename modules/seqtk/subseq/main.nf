process SEQTK_SUBSEQ {
    tag '$sequences'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
    } else {
        container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    }

    input:
    path sequences
    path filter_list

    output:
    path "*.gz"         , emit: sequences
    path "versions.yml" , emit: versions

    script:
    def prefix   = options.suffix ?: ''
    def ext = "fa"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fq"
    }
    """
    seqtk \\
        subseq \\
        $args \\
        $sequences \\
        $filter_list | \\
        gzip --no-name > ${sequences}${prefix}.${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
