process SICKLE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::sickle-trim=1.33"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sickle-trim:1.33--h7132678_7':
        'biocontainers/sickle-trim:1.33--h5bf99c6_6' }"

    input:
    tuple val(meta), path(reads), val(qual_type)

    output:
    tuple val(meta), path("${prefix}.se.trimmed.fastq.gz"),        optional:true, emit: single_trimmed
    tuple val(meta), path("${prefix}.pe{1,2}.trimmed.fastq.gz"),   optional:true, emit: paired_trimmed
    tuple val(meta), path("${prefix}.singleton.trimmed.fastq.gz"), optional:true, emit: singleton_trimmed
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sickle: \$(sickle --version|awk 'NR==1{print \$3}')
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sickle: \$(sickle --version|awk 'NR==1{print \$3}')
    END_VERSIONS
    """
    }
}
