process QCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qcat=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qcat:1.1.0--py_0' :
        'biocontainers/qcat:1.1.0--py_0' }"

    input:
    tuple val(meta), path(reads)
    val   barcode_kit

    output:
    tuple val(meta), path("fastq/*.fastq.gz"), emit: reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## Unzip fastq file
    ## qcat doesn't support zipped files yet
    FILE=$reads
    if [[ \$FILE == *.gz ]]
    then
        zcat $reads > unzipped.fastq
        FILE=unzipped.fastq
    fi

    qcat \\
        -f \$FILE \\
        -b ./fastq \\
        --kit $barcode_kit

    ## Zip fastq files
    gzip fastq/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qcat: \$(qcat --version 2>&1 | sed 's/^.*qcat //; s/ .*\$//')
    END_VERSIONS
    """
}
