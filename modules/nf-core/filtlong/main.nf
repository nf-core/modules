process FILTLONG {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::filtlong=0.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/filtlong:0.2.1--h9a82719_0' :
        'biocontainers/filtlong:0.2.1--h9a82719_0' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def short_reads = !shortreads ? "" : meta.single_end ? "-1 $shortreads" : "-1 ${shortreads[0]} -2 ${shortreads[1]}"
    if ("$longreads" == "${prefix}.fastq.gz") error "Longread FASTQ input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    filtlong \\
        $short_reads \\
        $args \\
        $longreads \\
        2> ${prefix}.log \\
        | gzip -n > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filtlong: \$( filtlong --version | sed -e "s/Filtlong v//g" )
    END_VERSIONS
    """
}
