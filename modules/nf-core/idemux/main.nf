process IDEMUX {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::idemux=0.1.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idemux:0.1.6--pyhdfd78af_0':
        'biocontainers/idemux:0.1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(samplesheet)

    output:
    tuple val(meta), path("[!undetermined]*.fastq.gz")          , emit: fastq
    tuple val(meta), path("undetermined_R?.fastq.gz")  , optional:true, emit: undetermined
    path "demultipexing_stats.tsv" , emit: stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    idemux \\
        --r1 ${reads[0]} \\
        --r2 ${reads[1]} \\
        --sample-sheet ${samplesheet} \\
        --out . \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        idemux: \$(idemux --version |& sed '1!d ; s/idemux //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    echo -e "sample_name\twritten_reads" > demultipexing_stats.tsv

    sed 1d ${samplesheet} | while IFS=, read -r sampleName _ _ _; do
        touch "\${sampleName}_R1.fastq"
        touch "\${sampleName}_R2.fastq"
        echo -e "\${sampleName}\t100" >> demultipexing_stats.tsv
    done

    touch undetermined_R1.fastq
    touch undetermined_R2.fastq
    gzip *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        idemux: \$(idemux --version |& sed '1!d ; s/idemux //')
    END_VERSIONS
    """
}
