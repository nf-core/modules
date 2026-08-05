process SICKLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sickle-trim:1.33--h7132678_7':
        'quay.io/biocontainers/sickle-trim:1.33--h5bf99c6_6' }"

    input:
    tuple val(meta), path(reads), val(qual_type)

    output:
    tuple val(meta), path("*.se.trimmed.fastq.gz"),        optional:true, emit: single_trimmed
    tuple val(meta), path("*.pe{1,2}.trimmed.fastq.gz"),   optional:true, emit: paired_trimmed
    tuple val(meta), path("*.singleton.trimmed.fastq.gz"), optional:true, emit: singleton_trimmed
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('sickle'), eval('sickle --version 2>&1 | head -1 | sed "s/sickle version //"'), topic: versions, emit: versions_sickle

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end){
    """
    sickle \\
        se \\
        $args \\
        -f $reads \\
        -t $qual_type \\
        -o ${prefix}.se.trimmed.fastq.gz \\
        -g \\
    >${prefix}.sickle.log
    """
    }
    else{
    """
    sickle \\
        pe \\
        $args \\
        -f ${reads[0]} \\
        -r ${reads[1]} \\
        -t $qual_type \\
        -o ${prefix}.pe1.trimmed.fastq.gz \\
        -p ${prefix}.pe2.trimmed.fastq.gz \\
        -s ${prefix}.singleton.trimmed.fastq.gz \\
        -g \\
    >${prefix}.sickle.log
    """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end){
    """
    echo "" | gzip > ${prefix}.se.trimmed.fastq.gz
    touch ${prefix}.sickle.log
    """
    }
    else {
    """
    echo "" | gzip > ${prefix}.pe1.trimmed.fastq.gz
    echo "" | gzip > ${prefix}.pe2.trimmed.fastq.gz
    echo "" | gzip > ${prefix}.singleton.trimmed.fastq.gz
    touch ${prefix}.sickle.log
    """
    }
}
