process FILTLONG {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::filtlong=0.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/filtlong:0.2.1--h9a82719_0' :
        'quay.io/biocontainers/filtlong:0.2.1--h9a82719_0' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path("${meta.id}_lr_filtlong.fastq.gz"), emit: reads
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def short_reads = meta.single_end ? "-1 $shortreads" : "-1 ${shortreads[0]} -2 ${shortreads[1]}"
    """
    filtlong \\
        $short_reads \\
        $args \\
        $longreads \\
        | gzip -n > ${prefix}_lr_filtlong.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filtlong: \$( filtlong --version | sed -e "s/Filtlong v//g" )
    END_VERSIONS
    """
}
