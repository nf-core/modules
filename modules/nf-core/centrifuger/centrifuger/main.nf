process CENTRIFUGER_CENTRIFUGER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuger:1.1.0--hf426362_0':
        'quay.io/biocontainers/centrifuger:1.1.0--hf426362_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(db)
    val save_unclassified
    val save_classified
    path barcode
    path umi

    output:
    tuple val(meta), path("*.tsv")                , emit: classification_file
    tuple val(meta), path("*.classified*.fq.gz")  , emit: fastq_classified  , optional: true
    tuple val(meta), path("*.unclassified*.fq.gz"), emit: fastq_unclassified, optional: true
    tuple val("${task.process}"), val('centrifuger'), eval("centrifuger -v 2>&1 | sed 's/Centrifuger v//'"),emit: versions_centrifuger,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-u ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    // Optional outputs
    def unclassified_arg = save_unclassified ? "--un ${prefix}.unclassified" : ""
    def classified_arg = save_classified ? "--cl ${prefix}.classified" : ""
    def barcode_arg = barcode ? "--barcode ${barcode}" : ""
    def umi_arg = umi ? "--UMI ${umi}" : ""

    """
    db_name=`find -L . -name "*.1.cfr" -not -name "._*"  | sed 's/\\.1.cfr\$//'`

    centrifuger \\
        -x \$db_name \\
        ${paired} \\
        ${unclassified_arg} \\
        ${classified_arg} \\
        ${barcode_arg} \\
        ${umi_arg} \\
        -t ${task.cpus} \\
        ${args} > ${prefix}.tsv
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    ## main output
    touch ${prefix}.tsv

    ## Optional outputs
    if ${save_unclassified}; then
        if ${meta.single_end}; then
            echo "" | gzip >  ${prefix}.unclassified.fq.gz
        else
            echo "" | gzip > ${prefix}.unclassified_1.fq.gz
            echo "" | gzip > ${prefix}.unclassified_2.fq.gz
        fi
    fi

    if ${save_classified}; then
        if ${meta.single_end}; then
            echo "" | gzip  > ${prefix}.classified.fq.gz
        else
           echo "" | gzip > ${prefix}.classified_1.fq.gz
           echo "" | gzip > ${prefix}.classified_2.fq.gz
        fi
    fi
    """
}
