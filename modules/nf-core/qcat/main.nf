process QCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qcat:1.1.0--py_0' :
        'quay.io/biocontainers/qcat:1.1.0--py_0' }"

    input:
    tuple val(meta), path(reads)
    val barcode_kit

    output:
    tuple val(meta), path("fastq/*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('qcat'), eval("qcat --version | sed 's/.*qcat //;s/ .*//'"), emit: versions_qcat, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Unzip fastq file
    ## qcat doesn't support zipped files yet
    FILE=${reads}
    if [[ \$FILE == *.gz ]]
    then
        zcat ${reads} > unzipped.fastq
        FILE=unzipped.fastq
    fi

    qcat \\
        -f \$FILE \\
        -b ./fastq \\
        --kit ${barcode_kit}

    ## Zip fastq files
    gzip fastq/*
    """

    stub:
    """
    mkdir -p fastq

    echo "" | gzip > fastq/barcode00.fastq.gz
    echo "" | gzip > fastq/none.fastq.gz
    """
}
